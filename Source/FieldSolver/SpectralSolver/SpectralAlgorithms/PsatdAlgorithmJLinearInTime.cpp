/* Copyright 2019
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PsatdAlgorithmJLinearInTime.H"

#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpX_Complex.H"

#include <AMReX_Array4.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BaseFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_PODVector.H>

#include <cmath>

#if WARPX_USE_PSATD

using namespace amrex::literals;

PsatdAlgorithmJLinearInTime::PsatdAlgorithmJLinearInTime(
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const SpectralFieldIndex& spectral_index,
    const int norder_x,
    const int norder_y,
    const int norder_z,
    const bool nodal,
    const amrex::Real dt,
    const bool time_averaging,
    const bool dive_cleaning,
    const bool divb_cleaning)
    // Initializer list
    : SpectralBaseAlgorithm(spectral_kspace, dm, spectral_index, norder_x, norder_y, norder_z, nodal),
    m_spectral_index(spectral_index),
    m_dt(dt),
    m_time_averaging(time_averaging),
    m_dive_cleaning(dive_cleaning),
    m_divb_cleaning(divb_cleaning)
{
    const amrex::BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Always allocate these coefficients
    C_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X1_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X2_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X3_coef = SpectralRealCoefficients(ba, dm, 1, 0);

    InitializeSpectralCoefficients(spectral_kspace, dm, dt);

    // Allocate these coefficients only with time averaging
    if (time_averaging)
    {
        X5_coef = SpectralRealCoefficients(ba, dm, 1, 0);
        X6_coef = SpectralRealCoefficients(ba, dm, 1, 0);
        InitializeSpectralCoefficientsAveraging(spectral_kspace, dm, dt);
    }
}

void
PsatdAlgorithmJLinearInTime::pushSpectralFields (SpectralFieldData& f) const
{
    const bool time_averaging = m_time_averaging;
    const bool dive_cleaning = m_dive_cleaning;
    const bool divb_cleaning = m_divb_cleaning;

    const amrex::Real dt = m_dt;
    const amrex::Real dt2 = dt*dt;

    const SpectralFieldIndex& Idx = m_spectral_index;

    // Loop over boxes
    for (amrex::MFIter mfi(f.fields); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> fields = f.fields[mfi].array();

        // These coefficients are always allocated
        amrex::Array4<const amrex::Real> C_arr = C_coef[mfi].array();
        amrex::Array4<const amrex::Real> S_ck_arr = S_ck_coef[mfi].array();
        amrex::Array4<const amrex::Real> X1_arr = X1_coef[mfi].array();
        amrex::Array4<const amrex::Real> X2_arr = X2_coef[mfi].array();
        amrex::Array4<const amrex::Real> X3_arr = X3_coef[mfi].array();

        amrex::Array4<const amrex::Real> X5_arr;
        amrex::Array4<const amrex::Real> X6_arr;
        if (time_averaging)
        {
            X5_arr = X5_coef[mfi].array();
            X6_arr = X6_coef[mfi].array();
        }

        // Extract pointers for the k vectors
        const amrex::Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if defined(WARPX_DIM_3D)
        const amrex::Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Record old values of the fields to be updated
            const Complex Ex_old = fields(i,j,k,Idx.Ex);
            const Complex Ey_old = fields(i,j,k,Idx.Ey);
            const Complex Ez_old = fields(i,j,k,Idx.Ez);
            const Complex Bx_old = fields(i,j,k,Idx.Bx);
            const Complex By_old = fields(i,j,k,Idx.By);
            const Complex Bz_old = fields(i,j,k,Idx.Bz);

            // Shortcuts for the values of J and rho
            const Complex Jx_old = fields(i,j,k,Idx.Jx);
            const Complex Jy_old = fields(i,j,k,Idx.Jy);
            const Complex Jz_old = fields(i,j,k,Idx.Jz);
            const Complex Jx_new = fields(i,j,k,Idx.Jx_new);
            const Complex Jy_new = fields(i,j,k,Idx.Jy_new);
            const Complex Jz_new = fields(i,j,k,Idx.Jz_new);
            const Complex rho_old = fields(i,j,k,Idx.rho_old);
            const Complex rho_mid = fields(i,j,k,Idx.rho_mid);
            const Complex rho_new = fields(i,j,k,Idx.rho_new);

            Complex F_old, G_old;
            if (dive_cleaning) F_old = fields(i,j,k,Idx.F);
            if (divb_cleaning) G_old = fields(i,j,k,Idx.G);

            // k vector values
            const amrex::Real kx = modified_kx_arr[i];
#if defined(WARPX_DIM_3D)
            const amrex::Real ky = modified_ky_arr[j];
            const amrex::Real kz = modified_kz_arr[k];
#else
            constexpr amrex::Real ky = 0._rt;
            const     amrex::Real kz = modified_kz_arr[j];
#endif
            // Physical constants and imaginary unit
            constexpr amrex::Real c2 = PhysConst::c * PhysConst::c;
            constexpr amrex::Real mu0 = PhysConst::mu0;
            constexpr amrex::Real ep0 = PhysConst::ep0;
            constexpr amrex::Real inv_ep0 = 1._rt / PhysConst::ep0;
            constexpr Complex I = Complex{0._rt, 1._rt};

            // These coefficients are initialized in the function InitializeSpectralCoefficients
            const amrex::Real C = C_arr(i,j,k);
            const amrex::Real S_ck = S_ck_arr(i,j,k);
            const amrex::Real X1 = X1_arr(i,j,k);
            const amrex::Real X2 = X2_arr(i,j,k);
            const amrex::Real X3 = X3_arr(i,j,k);
            const amrex::Real X4 = - S_ck / PhysConst::ep0;

            const amrex::Real kx2 = kx*kx;
            const amrex::Real ky2 = ky*ky;
            const amrex::Real kz2 = kz*kz;
            const amrex::Real knorm2 = kx2 + ky2 + kz2;
            const amrex::Real knorm4 = knorm2*knorm2;

            const Complex Jx_c1 = (Jx_new-Jx_old) / dt;
            const Complex Jy_c1 = (Jy_new-Jy_old) / dt;
            const Complex Jz_c1 = (Jz_new-Jz_old) / dt;

            if (!dive_cleaning && !divb_cleaning)
            {
                if (knorm2 == 0._rt)
                {
                    fields(i,j,k,Idx.Ex) = Ex_old - mu0*c2*dt*Jx_old - 0.5_rt*mu0*c2*dt2*Jx_c1;
                    fields(i,j,k,Idx.Ey) = Ey_old - mu0*c2*dt*Jy_old - 0.5_rt*mu0*c2*dt2*Jy_c1;
                    fields(i,j,k,Idx.Ez) = Ez_old - mu0*c2*dt*Jz_old - 0.5_rt*mu0*c2*dt2*Jz_c1;
                }
                else // knorm2 != 0
                {
                    const amrex::Real C1 = (kx2 + ky2*C + kz2*C) / knorm2;
                    const amrex::Real C2 = (ky2 + kz2*C + kx2*C) / knorm2;
                    const amrex::Real C3 = (kz2 + kx2*C + ky2*C) / knorm2;

                    const amrex::Real C4 = kx*ky*(1._rt-C) / knorm2;
                    const amrex::Real C5 = ky*kz*(1._rt-C) / knorm2;
                    const amrex::Real C6 = kz*kx*(1._rt-C) / knorm2;

                    const amrex::Real C7 = c2*kx*S_ck;
                    const amrex::Real C8 = c2*ky*S_ck;
                    const amrex::Real C9 = c2*kz*S_ck;

                    const amrex::Real C10 = -mu0*c2*(dt*kx2 + ky2*S_ck + kz2*S_ck) / knorm2;
                    const amrex::Real C11 = -mu0*c2*(dt*ky2 + kz2*S_ck + kx2*S_ck) / knorm2;
                    const amrex::Real C12 = -mu0*c2*(dt*kz2 + kx2*S_ck + ky2*S_ck) / knorm2;

                    const amrex::Real C13 = mu0*c2*kx*ky*(S_ck-dt) / knorm2;
                    const amrex::Real C14 = mu0*c2*ky*kz*(S_ck-dt) / knorm2;
                    const amrex::Real C15 = mu0*c2*kz*kx*(S_ck-dt) / knorm2;

                    const amrex::Real C16 = mu0*kx*(1._rt-C) / knorm2;
                    const amrex::Real C17 = mu0*ky*(1._rt-C) / knorm2;
                    const amrex::Real C18 = mu0*kz*(1._rt-C) / knorm2;

                    const amrex::Real C19 = -mu0*(c2*dt2*kx2*knorm2 + 2._rt*ky2*(1._rt-C)
                                            + 2._rt*kz2*(1._rt-C)) / (2._rt*knorm4);
                    const amrex::Real C20 = -mu0*(c2*dt2*ky2*knorm2 + 2._rt*kz2*(1._rt-C)
                                            + 2._rt*kx2*(1._rt-C)) / (2._rt*knorm4);
                    const amrex::Real C21 = -mu0*(c2*dt2*kz2*knorm2 + 2._rt*kx2*(1._rt-C)
                                            + 2._rt*ky2*(1._rt-C)) / (2._rt*knorm4);

                    const amrex::Real C22 = -mu0*kx*ky*(c2*dt2*knorm2
                                            + 2._rt*(C-1._rt)) / (2._rt*knorm4);
                    const amrex::Real C23 = -mu0*ky*kz*(c2*dt2*knorm2
                                            + 2._rt*(C-1._rt)) / (2._rt*knorm4);
                    const amrex::Real C24 = -mu0*kz*kx*(c2*dt2*knorm2
                                            + 2._rt*(C-1._rt)) / (2._rt*knorm4);

                    const amrex::Real C25 = mu0*kx*(S_ck-dt) / knorm2;
                    const amrex::Real C26 = mu0*ky*(S_ck-dt) / knorm2;
                    const amrex::Real C27 = mu0*kz*(S_ck-dt) / knorm2;

                    fields(i,j,k,Idx.Ex) = C1*Ex_old + C4*Ey_old + C6*Ez_old
                                           - I*C9*By_old + I*C8*Bz_old
                                           + C10*Jx_old + C13*Jy_old + C15*Jz_old;
                                           + C19*Jx_c1 + C22*Jy_c1 + C24*Jz_c1;

                    fields(i,j,k,Idx.Ey) = C4*Ex_old + C2*Ey_old + C5*Ez_old
                                           + I*C9*Bx_old - I*C7*Bz_old
                                           + C13*Jx_old + C11*Jy_old + C14*Jz_old
                                           + C22*Jx_c1 + C20*Jy_c1 + C23*Jz_c1;

                    fields(i,j,k,Idx.Ez) = C6*Ex_old + C5*Ey_old + C3*Ez_old
                                           - I*C8*Bx_old + I*C7*By_old
                                           + C15*Jx_old + C14*Jy_old + C12*Jz_old
                                           + C24*Jx_c1 + C23*Jy_c1 + C21*Jz_c1;

                    fields(i,j,k,Idx.Bx) = C1*Bx_old + C4*By_old + C6*Bz_old
                                           + I*C9/c2*Ey_old - I*C8/c2*Ez_old
                                           - I*C18*Jy_old + I*C17*Jz_old
                                           + I*C27*Jy_c1 - I*C26*Jz_c1;

                    fields(i,j,k,Idx.By) = C4*Bx_old + C2*By_old + C5*Bz_old
                                           - I*C9/c2*Ex_old + I*C7/c2*Ez_old
                                           + I*C18*Jx_old - I*C16*Jz_old
                                           - I*C27*Jx_c1 + I*C25*Jz_c1;

                    fields(i,j,k,Idx.Bz) = C6*Bx_old + C5*By_old + C3*Bz_old
                                           + I*C8/c2*Ex_old - I*C7/c2*Ey_old
                                           - I*C17*Jx_old + I*C16*Jy_old
                                           + I*C26*Jx_c1 - I*C25*Jy_c1;
                }
            }
            else if (dive_cleaning && divb_cleaning)
            {
                const Complex rho_c2 = (-4._rt*rho_mid + 2._rt*rho_new + 2._rt*rho_old) / dt2;
                const Complex rho_c1 = (4._rt*rho_mid - rho_new - 3._rt*rho_old) / dt;

                if (knorm2 == 0._rt)
                {
                    fields(i,j,k,Idx.Ex) = Ex_old - mu0*c2*dt*Jx_old - 0.5_rt*mu0*c2*dt2*Jx_c1;
                    fields(i,j,k,Idx.Ey) = Ey_old - mu0*c2*dt*Jy_old - 0.5_rt*mu0*c2*dt2*Jy_c1;
                    fields(i,j,k,Idx.Ez) = Ez_old - mu0*c2*dt*Jz_old - 0.5_rt*mu0*c2*dt2*Jz_c1;

                    fields(i,j,k,Idx.F) = F_old - 0.5_rt*mu0*c2*dt2*rho_c1 - inv_ep0*dt*rho_old;
                }
                else // knorm2 != 0
                {
                    const amrex::Real C1 = c2*kx*S_ck;
                    const amrex::Real C2 = c2*ky*S_ck;
                    const amrex::Real C3 = c2*kz*S_ck;

                    const amrex::Real C4 = mu0*c2*kx*(S_ck-dt) / knorm2;
                    const amrex::Real C5 = mu0*c2*ky*(S_ck-dt) / knorm2;
                    const amrex::Real C6 = mu0*c2*kz*(S_ck-dt) / knorm2;

                    const amrex::Real C7 = mu0*kx*(1._rt-C) / knorm2;
                    const amrex::Real C8 = mu0*ky*(1._rt-C) / knorm2;
                    const amrex::Real C9 = mu0*kz*(1._rt-C) / knorm2;

                    const amrex::Real C10 = -mu0*kx*(c2*dt2*knorm2
                                            + 2._rt*(C-1._rt)) / knorm4;
                    const amrex::Real C11 = -mu0*ky*(c2*dt2*knorm2
                                            + 2._rt*(C-1._rt)) / knorm4;
                    const amrex::Real C12 = -mu0*kz*(c2*dt2*knorm2
                                            + 2._rt*(C-1._rt)) / knorm4;

                    const amrex::Real C13 = mu0*kx*(S_ck-dt) / knorm2;
                    const amrex::Real C14 = mu0*ky*(S_ck-dt) / knorm2;
                    const amrex::Real C15 = mu0*kz*(S_ck-dt) / knorm2;

                    const amrex::Real C16 = mu0*(C-1._rt) / knorm2;
                    const amrex::Real C17 = 2._rt*mu0*(S_ck-dt) / knorm2;

                    fields(i,j,k,Idx.Ex) = C*Ex_old - I*C3*By_old + I*C2*Bz_old
                                           + I*C1*F_old - inv_ep0*S_ck*Jx_old + C16*Jx_c1
                                           + I*C10*rho_c2 + I*C4*rho_c1 - I*c2*C7*rho_old;

                    fields(i,j,k,Idx.Ey) = C*Ey_old + I*C3*Bx_old - I*C1*Bz_old
                                           + I*C2*F_old - inv_ep0*S_ck*Jy_old + C16*Jy_c1
                                           + I*C11*rho_c2 + I*C5*rho_c1 - I*c2*C8*rho_old;

                    fields(i,j,k,Idx.Ez) = C*Ez_old - I*C2*Bx_old + I*C1*By_old
                                           + I*C3*F_old - inv_ep0*S_ck*Jz_old + C16*Jz_c1
                                           + I*C12*rho_c2 + I*C6*rho_c1 - I*c2*C9*rho_old;

                    fields(i,j,k,Idx.Bx) = C*Bx_old + I*C3/c2*Ey_old - I*C2/c2*Ez_old
                                           + I*C1*G_old - I*C9*Jy_old + I*C8*Jz_old
                                           + I*C15*Jy_c1 - I*C14*Jz_c1;

                    fields(i,j,k,Idx.By) = C*By_old - I*C3/c2*Ex_old + I*C1/c2*Ez_old
                                           + I*C2*G_old + I*C9*Jx_old - I*C7*Jz_old
                                           - I*C15*Jx_c1 + I*C13*Jz_c1;

                    fields(i,j,k,Idx.Bz) = C*Bz_old + I*C2/c2*Ex_old - I*C1/c2*Ey_old
                                           + I*C3*G_old - I*C8*Jx_old + I*C7*Jy_old
                                           + I*C14*Jx_c1 - I*C13*Jy_c1;

                    fields(i,j,k,Idx.F) = C*F_old
                                          + I*C1/c2*Ex_old + I*C2/c2*Ey_old + I*C3/c2*Ez_old
                                          - I*C7*Jx_old - I*C8*Jy_old - I*C9*Jz_old
                                          + I*C13*Jx_c1 + I*C14*Jy_c1 + I*C15*Jz_c1
                                          + C17*rho_c2 + C16*rho_c1 - inv_ep0*S_ck*rho_old;

                    fields(i,j,k,Idx.G) = C*G_old
                                          + I*C1/c2*Bx_old + I*C2/c2*By_old + I*C3/c2*Bz_old;
                }
            }

            //// Update equations for E in the formulation with rho

            //fields(i,j,k,Idx.Ex) = C * Ex_old
            //    + I * c2 * S_ck * (ky * Bz_old - kz * By_old)
            //    + X4 * Jx_old - I * (X2 * rho_new - X3 * rho_old) * kx - X1 * (Jx_new - Jx_old) / dt;

            //fields(i,j,k,Idx.Ey) = C * Ey_old
            //    + I * c2 * S_ck * (kz * Bx_old - kx * Bz_old)
            //    + X4 * Jy_old - I * (X2 * rho_new - X3 * rho_old) * ky - X1 * (Jy_new - Jy_old) / dt;

            //fields(i,j,k,Idx.Ez) = C * Ez_old
            //    + I * c2 * S_ck * (kx * By_old - ky * Bx_old)
            //    + X4 * Jz_old - I * (X2 * rho_new - X3 * rho_old) * kz - X1 * (Jz_new - Jz_old) / dt;

            //// Update equations for B

            //fields(i,j,k,Idx.Bx) = C * Bx_old
            //    - I * S_ck * (ky * Ez_old - kz * Ey_old) + I * X1 * (ky * Jz_old - kz * Jy_old)
            //    + I * X2/c2 * (ky * (Jz_new - Jz_old) - kz * (Jy_new - Jy_old));

            //fields(i,j,k,Idx.By) = C * By_old
            //    - I * S_ck * (kz * Ex_old - kx * Ez_old) + I * X1 * (kz * Jx_old - kx * Jz_old)
            //    + I * X2/c2 * (kz * (Jx_new - Jx_old) - kx * (Jz_new - Jz_old));

            //fields(i,j,k,Idx.Bz) = C * Bz_old
            //    - I * S_ck * (kx * Ey_old - ky * Ex_old) + I * X1 * (kx * Jy_old - ky * Jx_old)
            //    + I * X2/c2 * (kx * (Jy_new - Jy_old) - ky * (Jx_new - Jx_old));

            //if (dive_cleaning)
            //{
            //    const Complex k_dot_E = kx * Ex_old + ky * Ey_old + kz * Ez_old;
            //    const Complex k_dot_J  = kx * Jx_old + ky * Jy_old + kz * Jz_old;
            //    const Complex k_dot_dJ = kx * (Jx_new - Jx_old) + ky * (Jy_new - Jy_old) + kz * (Jz_new - Jz_old);

            //    fields(i,j,k,Idx.Ex) += I * c2 * S_ck * F_old * kx;
            //    fields(i,j,k,Idx.Ey) += I * c2 * S_ck * F_old * ky;
            //    fields(i,j,k,Idx.Ez) += I * c2 * S_ck * F_old * kz;

            //    fields(i,j,k,Idx.F) = C * F_old + S_ck * (I * k_dot_E - rho_old * inv_ep0)
            //        - X1 * ((rho_new - rho_old) / dt + I * k_dot_J) - I * X2/c2 * k_dot_dJ;
            //}

            //if (divb_cleaning)
            //{
            //    const Complex k_dot_B = kx * Bx_old + ky * By_old + kz * Bz_old;

            //    fields(i,j,k,Idx.Bx) += I * S_ck * G_old * kx;
            //    fields(i,j,k,Idx.By) += I * S_ck * G_old * ky;
            //    fields(i,j,k,Idx.Bz) += I * S_ck * G_old * kz;

            //    fields(i,j,k,Idx.G) = C * G_old + I * c2 * S_ck * k_dot_B;
            //}

            if (time_averaging)
            {
                const amrex::Real X5 = X5_arr(i,j,k);
                const amrex::Real X6 = X6_arr(i,j,k);

                // TODO: Here the code is *accumulating* the average,
                // because it is meant to be used with sub-cycling
                // maybe this should be made more generic

                fields(i,j,k,Idx.Ex_avg) += S_ck * Ex_old
                    + I * c2 * ep0 * X1 * (ky * Bz_old - kz * By_old)
                    + I * X5 * rho_old * kx + I * X6 * rho_new * kx + X3/c2 * Jx_old - X2/c2 * Jx_new;

                fields(i,j,k,Idx.Ey_avg) += S_ck * Ey_old
                    + I * c2 * ep0 * X1 * (kz * Bx_old - kx * Bz_old)
                    + I * X5 * rho_old * ky + I * X6 * rho_new * ky + X3/c2 * Jy_old - X2/c2 * Jy_new;

                fields(i,j,k,Idx.Ez_avg) += S_ck * Ez_old
                    + I * c2 * ep0 * X1 * (kx * By_old - ky * Bx_old)
                    + I * X5 * rho_old * kz + I * X6 * rho_new * kz + X3/c2 * Jz_old - X2/c2 * Jz_new;

                fields(i,j,k,Idx.Bx_avg) += S_ck * Bx_old
                    - I * ep0 * X1 * (ky * Ez_old - kz * Ey_old)
                    - I * X5/c2 * (ky * Jz_old - kz * Jy_old) - I * X6/c2 * (ky * Jz_new - kz * Jy_new);

                fields(i,j,k,Idx.By_avg) += S_ck * By_old
                    - I * ep0 * X1 * (kz * Ex_old - kx * Ez_old)
                    - I * X5/c2 * (kz * Jx_old - kx * Jz_old) - I * X6/c2 * (kz * Jx_new - kx * Jz_new);

                fields(i,j,k,Idx.Bz_avg) += S_ck * Bz_old
                    - I * ep0 * X1 * (kx * Ey_old - ky * Ex_old)
                    - I * X5/c2 * (kx * Jy_old - ky * Jx_old) - I * X6/c2 * (kx * Jy_new - ky * Jx_new);

                if (dive_cleaning)
                {
                    fields(i,j,k,Idx.Ex_avg) += I * c2 * ep0 * X1 * F_old * kx;
                    fields(i,j,k,Idx.Ey_avg) += I * c2 * ep0 * X1 * F_old * ky;
                    fields(i,j,k,Idx.Ez_avg) += I * c2 * ep0 * X1 * F_old * kz;
                }

                if (divb_cleaning)
                {
                    fields(i,j,k,Idx.Bx_avg) += I * ep0 * X1 * G_old * kx;
                    fields(i,j,k,Idx.By_avg) += I * ep0 * X1 * G_old * ky;
                    fields(i,j,k,Idx.Bz_avg) += I * ep0 * X1 * G_old * kz;
                }
            }
        });
    }
}

void PsatdAlgorithmJLinearInTime::InitializeSpectralCoefficients (
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const amrex::Real dt)
{
    const amrex::BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Loop over boxes and allocate the corresponding coefficients for each box
    for (amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const amrex::Real* kx_s = modified_kx_vec[mfi].dataPtr();
#if defined(WARPX_DIM_3D)
        const amrex::Real* ky_s = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* kz_s = modified_kz_vec[mfi].dataPtr();

        // Coefficients always allocated
        amrex::Array4<amrex::Real> C = C_coef[mfi].array();
        amrex::Array4<amrex::Real> S_ck = S_ck_coef[mfi].array();
        amrex::Array4<amrex::Real> X1 = X1_coef[mfi].array();
        amrex::Array4<amrex::Real> X2 = X2_coef[mfi].array();
        amrex::Array4<amrex::Real> X3 = X3_coef[mfi].array();

        // Loop over indices within one box
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of k vector
            const amrex::Real knorm_s = std::sqrt(
                std::pow(kx_s[i], 2) +
#if defined(WARPX_DIM_3D)
                std::pow(ky_s[j], 2) + std::pow(kz_s[k], 2));
#else
                std::pow(kz_s[j], 2));
#endif
            // Physical constants and imaginary unit
            constexpr amrex::Real c = PhysConst::c;
            constexpr amrex::Real ep0 = PhysConst::ep0;

            const amrex::Real c2 = std::pow(c, 2);
            const amrex::Real dt2 = std::pow(dt, 2);

            const amrex::Real om_s = c * knorm_s;
            const amrex::Real om2_s = std::pow(om_s, 2);

            // C
            C(i,j,k) = std::cos(om_s * dt);

            // S_ck
            if (om_s != 0.)
            {
                S_ck(i,j,k) = std::sin(om_s * dt) / om_s;
            }
            else // om_s = 0
            {
                S_ck(i,j,k) = dt;
            }

            // X1 (multiplies i*([k] \times J) in the update equation for update B)
            if (om_s != 0.)
            {
                X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0 * om2_s);
            }
            else // om_s = 0
            {
                X1(i,j,k) = 0.5_rt * dt2 / ep0;
            }

            // X2 (multiplies rho_new in the update equation for E)
            if (om_s != 0.)
            {
                X2(i,j,k) = c2 * (dt - S_ck(i,j,k)) / (ep0 * dt * om2_s);
            }
            else // om_s = 0
            {
                X2(i,j,k) = c2 * dt2 / (6._rt * ep0);
            }

            // X3 (multiplies rho_old in the update equation for E)
            if (om_s != 0.)
            {
                X3(i,j,k) = c2 * (dt * C(i,j,k) - S_ck(i,j,k)) / (ep0 * dt * om2_s);
            }
            else // om_s = 0
            {
                X3(i,j,k) = - c2 * dt2 / (3._rt * ep0);
            }
        });
    }
}

void PsatdAlgorithmJLinearInTime::InitializeSpectralCoefficientsAveraging (
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const amrex::Real dt)
{
    const amrex::BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Loop over boxes and allocate the corresponding coefficients for each box
    for (amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const amrex::Real* kx_s = modified_kx_vec[mfi].dataPtr();
#if defined(WARPX_DIM_3D)
        const amrex::Real* ky_s = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* kz_s = modified_kz_vec[mfi].dataPtr();

        amrex::Array4<amrex::Real const> C = C_coef[mfi].array();
        amrex::Array4<amrex::Real const> S_ck = S_ck_coef[mfi].array();

        amrex::Array4<amrex::Real> X5 = X5_coef[mfi].array();
        amrex::Array4<amrex::Real> X6 = X6_coef[mfi].array();

        // Loop over indices within one box
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of k vector
            const amrex::Real knorm_s = std::sqrt(
                std::pow(kx_s[i], 2) +
#if defined(WARPX_DIM_3D)
                std::pow(ky_s[j], 2) + std::pow(kz_s[k], 2));
#else
                std::pow(kz_s[j], 2));
#endif
            // Physical constants and imaginary unit
            constexpr amrex::Real c = PhysConst::c;
            constexpr amrex::Real c2 = c*c;
            constexpr amrex::Real ep0 = PhysConst::ep0;

            // Auxiliary coefficients
            const amrex::Real dt3 = dt * dt * dt;

            const amrex::Real om_s  = c * knorm_s;
            const amrex::Real om2_s = om_s * om_s;
            const amrex::Real om4_s = om2_s * om2_s;

            if (om_s != 0.)
            {
                X5(i,j,k) = c2 / ep0 * (S_ck(i,j,k) / om2_s - (1._rt - C(i,j,k)) / (om4_s * dt)
                                        - 0.5_rt * dt / om2_s);
            }
            else
            {
                X5(i,j,k) = - c2 * dt3 / (8._rt * ep0);
            }

            if (om_s != 0.)
            {
                X6(i,j,k) = c2 / ep0 * ((1._rt - C(i,j,k)) / (om4_s * dt) - 0.5_rt * dt / om2_s);
            }
            else
            {
                X6(i,j,k) = - c2 * dt3 / (24._rt * ep0);
            }
        });
    }
}

void PsatdAlgorithmJLinearInTime::CurrentCorrection (SpectralFieldData& field_data)
{
    // Profiling
    BL_PROFILE("PsatdAlgorithmJLinearInTime::CurrentCorrection");

    amrex::ignore_unused(field_data);
    amrex::Abort(Utils::TextMsg::Err(
        "Current correction not implemented for multi-J PSATD algorithm"));
}

void
PsatdAlgorithmJLinearInTime::VayDeposition (SpectralFieldData& field_data)
{
    // Profiling
    BL_PROFILE("PsatdAlgorithmJLinearInTime::VayDeposition()");

    amrex::ignore_unused(field_data);
    amrex::Abort(Utils::TextMsg::Err(
        "Vay deposition not implemented for multi-J PSATD algorithm"));
}

#endif // WARPX_USE_PSATD

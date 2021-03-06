/* Copyright 2019 Axel Huebl, Maxence Thevenet, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <WarpX.H>

#include <AMReX_BaseFab.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_Config.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_SPACE.H>
#include <AMReX_TagBox.H>

#include <AMReX_BaseFwd.H>

using namespace amrex;

void
WarpX::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    const Real* problo = Geom(lev).ProbLo();
    const Real* dx = Geom(lev).CellSize();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tags); mfi.isValid(); ++mfi)
    {
        auto& fab = tags[mfi];
        const Box& bx = fab.box();
        for (BoxIterator bi(bx); bi.ok(); ++bi)
        {
            const IntVect& cell = bi();
            RealVect pos {AMREX_D_DECL((cell[0]+0.5_rt)*dx[0]+problo[0],
                                       (cell[1]+0.5_rt)*dx[1]+problo[1],
                                       (cell[2]+0.5_rt)*dx[2]+problo[2])};
            if (pos > fine_tag_lo && pos < fine_tag_hi) {
                fab(cell) = TagBox::SET;
            }
        }
    }
}

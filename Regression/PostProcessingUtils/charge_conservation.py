# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import numpy as np
from scipy.constants import e, epsilon_0


def _get_input_bool(input_dict, key, expected_value):
    # Input parameters are stored as lists of parsed string values. Return True only
    # when the requested key is present and its first value matches the expected value.
    value = input_dict.get(key)
    return value is not None and value[0] == expected_value


def _get_field_array(data, field):
    # Analysis scripts pass yt fields from different frontends. Convert the field
    # to a plain NumPy array while preserving the old behavior of each script.
    field_data = data[field]
    if hasattr(field_data, "to_ndarray"):
        return field_data.to_ndarray()
    if hasattr(field_data, "value"):
        return field_data.value
    if hasattr(field_data, "v"):
        return field_data.v
    return np.asarray(field_data)


def _default_check_options(input_file):
    from input_file_parser import parse_input_file

    # These flags determine whether the charge conservation check should run
    # and whether tolerances need to be relaxed.
    input_dict = parse_input_file(input_file)
    geometry_dims_rz = _get_input_bool(input_dict, "geometry.dims", "RZ")
    current_correction = _get_input_bool(input_dict, "psatd.current_correction", "1")
    current_deposition_vay = _get_input_bool(
        input_dict, "algo.current_deposition", "vay"
    )
    current_deposition_esirkepov = _get_input_bool(
        input_dict, "algo.current_deposition", "esirkepov"
    )
    maxwell_solver_psatd = _get_input_bool(input_dict, "algo.maxwell_solver", "psatd")

    # Check with current correction, Vay current deposition, and Esirkepov current deposition.
    # Do not check Esirkepov in RZ geometry, since that combination currently produces larger
    # numerical errors that need further investigation. Also do not check Esirkepov combined
    # with PSATD, since that combination does not conserve charge except for spectral order 2.
    do_check = (
        (
            current_deposition_esirkepov
            and not (geometry_dims_rz or maxwell_solver_psatd)
        )
        or current_correction
        or current_deposition_vay
    )

    # Default tolerance for the infinity-norm of the relative error between div(E) and rho/eps0.
    # This is relaxed for certain deposition schemes that produce larger numerical error.
    tolerance = 1.0e-11
    if current_correction:
        tolerance = 1.0e-9
    elif current_deposition_vay:
        tolerance = 1.0e-3

    return do_check, tolerance


def check_charge_conservation(
    data,
    *,
    tolerance=None,
    do_check=None,
    norm="relative_linf",
    normalization=None,
    trim=None,
    rho_field=("boxlib", "rho"),
    divE_field=("boxlib", "divE"),
    input_file="./warpx_used_inputs",
    title="Check charge conservation:",
    error_label=None,
    tolerance_label="tolerance",
):
    """Check Gauss's law using plotfile fields."""
    # If no explicit settings are passed, use the historical Langmuir behavior:
    # parse warpx_used_inputs to decide whether to run and which tolerance to use.
    # Other analyses pass explicit tolerances and select the norm that reproduces
    # their previous charge conservation check.
    if do_check is None and tolerance is None:
        do_check, tolerance = _default_check_options(input_file)
    else:
        if do_check is None:
            do_check = True
        if tolerance is None:
            tolerance = 1.0e-11

    if not do_check:
        return None

    # Load rho and div(E) from the diagnostic output. The field names are configurable
    # so callers can reuse this function with different plotfile frontend names.
    rho = _get_field_array(data, rho_field)
    divE = _get_field_array(data, divE_field)

    # Compute the requested error norm. The relative norms compare div(E) to rho/eps0,
    # while the normalized norms measure the local Gauss-law residual normalized by e*n0.
    if norm == "relative_linf":
        rho_over_eps0 = rho / epsilon_0
        error = np.amax(np.abs(divE - rho_over_eps0)) / np.amax(np.abs(rho_over_eps0))
    elif norm == "relative_linf_max":
        rho_over_eps0 = rho / epsilon_0
        error = np.amax(np.abs(divE - rho_over_eps0)) / max(
            np.amax(divE), np.amax(rho_over_eps0)
        )
    elif norm in ("normalized_rms", "normalized_linf"):
        if normalization is None:
            raise ValueError("normalization is required for normalized charge norms")
        drho = (rho - epsilon_0 * divE) / e / normalization
        if trim is not None:
            drho = drho[trim]
        if norm == "normalized_rms":
            error = np.sqrt((drho**2).sum() / drho.size)
        else:
            error = np.abs(drho).max()
    else:
        raise ValueError(f"Unknown charge conservation norm: {norm}")

    # Use default output labels that match the original analysis scripts, unless the
    # caller provides custom labels for its printed diagnostics.
    if error_label is None:
        error_label = {
            "relative_linf": "error_rel",
            "relative_linf_max": "err_charge",
            "normalized_rms": "drho_rms",
            "normalized_linf": "drho_max",
        }[norm]

    print(title)
    print(f"{error_label} = {error}")
    print(f"{tolerance_label} = {tolerance}")
    assert error < tolerance

    return error

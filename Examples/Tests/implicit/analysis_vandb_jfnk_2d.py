#!/usr/bin/env python3

# Copyright 2024 Justin Angus
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from the script `inputs_vandb_2d`.
# This simulates a 2D periodic plasma using the implicit solver
# with the Villasenor deposition using shape factor 2.
import sys

import numpy as np
import yt
from charge_conservation import check_charge_conservation

field_energy = np.loadtxt("diags/reducedfiles/field_energy.txt", skiprows=1)
particle_energy = np.loadtxt("diags/reducedfiles/particle_energy.txt", skiprows=1)

total_energy = field_energy[:, 2] + particle_energy[:, 2]

delta_E = (total_energy - total_energy[0]) / total_energy[0]
max_delta_E = np.abs(delta_E).max()

# This case should have near machine precision conservation of energy
tolerance_rel_energy = 2.0e-14
tolerance_rel_charge = 2.0e-15

print(f"max change in energy: {max_delta_E}")
print(f"tolerance: {tolerance_rel_energy}")

assert max_delta_E < tolerance_rel_energy

# check for machine precision conservation of charge density
n0 = 1.0e30

pltdir = sys.argv[1]
ds = yt.load(pltdir)
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)

check_charge_conservation(
    data,
    tolerance=tolerance_rel_charge,
    norm="normalized_rms",
    normalization=n0,
    title="rms error in charge conservation:",
)

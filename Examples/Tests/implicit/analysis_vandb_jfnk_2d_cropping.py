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

import yt
from charge_conservation import check_charge_conservation

# check for machine precision conservation of charge density
n0 = 1.0e12

pltdir = sys.argv[1]
ds = yt.load(pltdir)
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)

tolerance_max_charge = 1.0e-13

check_charge_conservation(
    data,
    tolerance=tolerance_max_charge,
    norm="normalized_linf",
    normalization=n0,
    title="max error in charge conservation:",
)

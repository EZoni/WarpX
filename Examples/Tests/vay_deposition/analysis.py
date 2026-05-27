#!/usr/bin/env python3

# Copyright 2019-2022
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import sys

import yt
from charge_conservation import check_charge_conservation

yt.funcs.mylog.setLevel(50)

# Plotfile data set
fn = sys.argv[1]
ds = yt.load(fn)

# Check relative L-infinity spatial norm of rho/epsilon_0 - div(E)
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)
tolerance = 1e-3
check_charge_conservation(
    data,
    tolerance=tolerance,
    title="Error on charge conservation:",
)

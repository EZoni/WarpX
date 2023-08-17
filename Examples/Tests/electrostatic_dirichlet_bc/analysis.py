#!/usr/bin/env python3

# Copyright 2021 Roelof Groenewald

# This script tests the time dependent Dirichlet boundary
# conditions in a 2D electrostatic simulation. An empty
# domain of 64 x 8 cells is simulated with periodic boundary
# conditions in the x directions and Dirichlet boundary
# conditions in the y direction with specified potentials
# of sine waves with different periods on the lo and hi side.
# One period of the hi side sine wave is simulated and the
# potentials at the boundaries compared to expectation.

# Possible running time: ~ 19 s

import glob
import os
import sys

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI
import numpy as np
import yt

plotfiles = sorted(glob.glob('dirichletbc_plt*'))[1:]
if len(plotfiles) == 0:
    plotfiles = sorted(glob.glob('Python_dirichletbc_plt*'))[1:]
assert len(plotfiles) > 0
opmdfile = './diags/diag2'

times = np.ones(len(plotfiles))
potentials_lo = np.zeros(len(plotfiles))
potentials_hi = np.zeros(len(plotfiles))

for ii, plotfile in enumerate(plotfiles):
    ds = yt.load(plotfile)
    times[ii] = (
        ds.current_time.item()
    )
    data = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
    )
    potentials_lo[ii] = np.mean(data['phi'].to_ndarray()[0])
    potentials_hi[ii] = np.mean(data['phi'].to_ndarray()[-1])

expected_potentials_lo = 150.0 * np.sin(2.0 * np.pi * 6.78e6 * times)
expected_potentials_hi = 450.0 * np.sin(2.0 * np.pi * 13.56e6 * times)

assert np.allclose(potentials_lo, expected_potentials_lo, rtol=0.1)
assert np.allclose(potentials_hi, expected_potentials_hi, rtol=0.1)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, output_file=plotfile, output_format='plotfile')
checksumAPI.evaluate_checksum(test_name, output_file=opmdfile, output_format='openpmd')

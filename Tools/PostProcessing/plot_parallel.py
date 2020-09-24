# Copyright 2019 David Grote, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import os
import glob
import matplotlib
import sys
import argparse
import yt
yt.funcs.mylog.setLevel(50)
import numpy as np
import matplotlib.pyplot as plt

'''
This script loops over all WarpX plotfiles in a directory and, for each
plotfile, saves a plot showing fields and particles.

Requires Python3 and that the yt version is higher than 3.5

It can be run serial
> python plot_parallel.py --path <path/to/plt/files> --serial
or parallel
> mpiexec -np 32 python plot_parallel.py --path <path/to/plt/files>

When running parallel, the plotfiles are distributed as evenly as possible between MPI ranks.

This script also proposes options to plot quantities over all time steps.
The quantity from all plotfiles is gathered to rank 0, and the evolution is plotted.
The default operation is the max of the given field. For example:
- "--plot_evolution Ey" will plot the max of Ey
- "--plot_particle_evolution species" for a given species will plot the RMS x versus the average z

To get help, execute
> python plot_parallel --help
'''

# Parse command line for options
parser = argparse.ArgumentParser()
parser.add_argument('--path', default = None,
                    help = 'path to plotfiles, defaults to "diags/". Plotfiles names must be "diag?????"',
                    metavar = '')
parser.add_argument('--plots_dir', default=None,
                    help = 'path where plots are saved, defaults to "diags/" or path if specified',
                    metavar = '')
parser.add_argument('--plotlib', default = 'yt', choices = ['yt', 'matplotlib'],
                    help = 'plotting library to use, defaults to yt',
                    metavar = '')
parser.add_argument('--field', default = 'Ez',
                    help = 'which grid field to plot (the central slice in y is plotted)',
                    metavar = '')
parser.add_argument('--pjump', default = 20,
                    help = 'when plotlib = matplotlib, we plot every pjump particle',
                    metavar = '')
parser.add_argument('--vmax', type = float, default = None,
                    help = 'if specified, the colormap will have bounds [-vmax, vmax]',
                    metavar = '')
parser.add_argument('--slicewidth', default = 10.e-6,
                    help = 'only particles with -slicewidth/2 < y < slicewidth/2 are plotted',
                    metavar = '')
parser.add_argument('--serial', action = 'store_true', default = False,
                    help = 'specifies running in serial, avoiding the import of MPI')
parser.add_argument('--species', dest = 'pslist', nargs = '+', type = str, default = None,
                    help = 'species to be plotted. By default, all species in the simulation are shown',
                    metavar = '')
parser.add_argument('--plot_evolution', type = str, default = None,
                    help = 'quantity to plot the evolution of across all data files',
                    metavar = '')
parser.add_argument('--plot_particle_evolution', type = str, default = None,
                    help = 'will plot the RMS x versus average z of the particles in the given species',
                    metavar = '')
args = parser.parse_args()

path = args.path
plots_dir = args.plots_dir
plotlib = args.plotlib
vmax = args.vmax
plot_evolution = args.plot_evolution
plot_particle_evolution = args.plot_particle_evolution

if path is None:
    path = 'diags/'
if plots_dir is None:
    plots_dir = path

# Sanity check
if int(sys.version[0]) != 3:
    print('WARNING: parallel analysis was tested only with Python3')

#matplotlib.rcParams.update({'font.size': 14})
pscolor = ['r','g','b','k','m','c','y','w']
pssize = 1.
# For 2D data, plot 2D array.
# For 3D data, plot central x-z slice.
yt_slicedir = {2:2, 3:1}

# Get list of particle species.
def get_species(a_file_list):
    # if user-specified, just return the user list
    if args.pslist is not None:
        return args.pslist
    # otherwise, loop over all plotfiles to get particle species list
    psset = set()
    for filename in a_file_list:
        ds = yt.load(filename)
        psset.add(ps for ps in ds.particle_types if ps != 'all')
    pslist = list(psset)
    pslist.sort()
    return pslist

def plot_snapshot(filename):

    print(filename)

    # Load plotfile
    ds = yt.load(filename)

    # Get number of dimension
    dim = ds.dimensionality

    # Plot grid fields
    if plotlib == 'matplotlib':
        fg = plt.figure(dpi = 200)
        ax = fg.add_subplot()
        # Read field quantities from yt dataset
        all_data_level_0 = ds.covering_grid(level = 0, left_edge = ds.domain_left_edge, dims = ds.domain_dimensions)
        F = all_data_level_0['boxlib', args.field].v.squeeze()
        if dim == 3:
            # FIXME Check cell-centering
            F = F[:,int(F.shape[1]+.5)//2,:]
        extent = [ds.domain_left_edge[dim-1], ds.domain_right_edge[dim-1],
                  ds.domain_left_edge[0], ds.domain_right_edge[0]]
        # Plot field quantities with matplotlib
        im = ax.imshow(F, aspect = 'auto', extent = extent, origin = 'lower')
        #for cl in im.collections:
        #    cl.set_edgecolor('face')
        if vmax is not None:
            im.set_clim(-vmax, vmax)

    if plotlib == 'yt':
        # Directly plot with yt
        if dim == 2:
            aspect = ds.domain_width[0] / ds.domain_width[1]
        if dim == 3:
            aspect = ds.domain_width[2] / ds.domain_width[0]
        sl = yt.SlicePlot(ds, yt_slicedir[dim], args.field, aspect = aspect)
        if vmax is not None:
            sl.set_zlim(-vmax, vmax)

    # Plot particle quantities
    for ispecies, pspecies in enumerate(pslist):
        if pspecies in [x[0] for x in ds.field_list]:

            if plotlib == 'matplotlib':
                # Read particle quantities from yt dataset
                ad = ds.all_data()
                xp = ad[pspecies, 'particle_position_x'].v
                if dim == 3:
                    yp = ad[pspecies, 'particle_position_y'].v
                    zp = ad[pspecies, 'particle_position_z'].v
                    select = yp**2<(args.slicewidth/2)**2
                    xp = xp[select] ; yp = yp[select] ; zp = zp[select]
                if dim == 2:
                    zp = ad[pspecies, 'particle_position_y'].v
                # Select randomly one every pjump particles
                random_indices = np.random.choice(xp.shape[0], int(xp.shape[0]/args.pjump))
                if dim == 2:
                    xp=xp[random_indices] ; zp=zp[random_indices]
                if dim == 3:
                    xp=xp[random_indices] ; yp=yp[random_indices] ; zp=zp[random_indices]
                ax.scatter(zp, xp, c = pscolor[ispecies], s = pssize, linewidth = pssize, marker = ',')

            if plotlib == 'yt':
                # Directly plot particles with yt
                sl.annotate_particles(width = (args.slicewidth, 'm'), p_size = pssize,
                                      ptype = pspecies, col = pscolor[ispecies])
    # Add labels to plot and save
    iteration = int(filename[-5:])
    image_file_name = os.path.join(plots_dir, 'plt_{:s}_{:s}_{:05d}.png'.format(args.field, plotlib, iteration))

    if plotlib == 'matplotlib':
        ax.set_xlabel('z (m)')
        ax.set_ylabel('x (m)')
        ax.set_title('{:s} at iteration {:d}, t = {:.2e} s'.format(args.field, iteration, ds.current_time.to_ndarray().mean()))
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        cb = fg.colorbar(im, ax = ax)
        cb.ax.tick_params()
        cb.ax.yaxis.get_offset_text().set()
        cb.ax.yaxis.offsetText.set_ha('center')
        cb.ax.yaxis.offsetText.set_va('bottom')
        cb.formatter.set_powerlimits((0,0)) # this does not work in log scale
        cb.update_ticks()
        ax.set_rasterized(True)
        fg.tight_layout()
        fg.savefig(image_file_name, dpi = 200)
        plt.close()

    if plotlib == 'yt':
        sl.annotate_grids()
        sl.save(image_file_name)

# Compute the evolved quantity from plotfile filename
def get_evolution_quantity(filename, quantity_name):
    # Load plotfile
    ds = yt.load( filename )
    # Get number of dimension
    dim = ds.dimensionality
    # Read field quantities from yt dataset
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F = all_data_level_0['boxlib', quantity_name].v.squeeze()
    zwin = (ds.domain_left_edge[dim-1]+ds.domain_right_edge[dim-1])/2
    quantity = np.amax(F)
    return zwin, quantity

def plot_evolved_quantity(zwin_arr, maxF_arr):
    plt.figure()
    plt.plot(zwin_arr, maxF_arr)
    plt.xlabel('z (m)')
    plt.ylabel('%s (S.I.)'%plot_evolution)
    plt.title('Field max evolution')
    plt.savefig(os.path.join(plots_dir, 'max_%s_evolution.pdf'%plot_evolution), bbox_inches='tight')

# Compute the evolved particle quantity from plotfile filename
def get_particle_evolution_quantity(filename, species):
    # Load plotfile
    ds = yt.load( filename )
    dim = ds.dimensionality
    ad = ds.all_data()
    x = ad[species, 'particle_position_x']
    if dim == 2:
        z = ad[species, 'particle_position_y']
    else:
        z = ad[species, 'particle_position_z']
    return np.mean(z), np.std(x)

def plot_particle_evolved_quantity(zbar, xstd):
    plt.figure()
    plt.plot(zbar, xstd)
    plt.xlabel('ave z (m)')
    plt.ylabel('rms x (m)')
    plt.title('%s evolution'%plot_particle_evolution)
    plt.savefig(os.path.join(plots_dir, '%s_evolution.pdf'%plot_particle_evolution), bbox_inches='tight')

def reduce_evolved_quantity(z, q):
    if size > 1:
        global_z = np.empty_like(z)
        global_q = np.empty_like(q)
        comm_world.Reduce(z, global_z, op=MPI.MAX)
        comm_world.Reduce(q, global_q, op=MPI.MAX)
        return global_z, global_q
    else:
        return z, q

### Analysis ###

# Get list of plotfiles
file_list = glob.glob(os.path.join(path, 'diag?????'))
file_list.sort()
nfiles = len(file_list)

# Get list of particle speciess to plot
pslist = get_species(file_list);

rank = 0
size = 1
if not args.serial:
    try:
        from mpi4py import MPI
        comm_world = MPI.COMM_WORLD
        rank = comm_world.Get_rank()
        size = comm_world.Get_size()
    except ImportError:
        pass

if rank == 0:
    print('Number of MPI ranks: %d'%size)
    print('Number of plotfiles: %s'%nfiles)
    print('List of species: ', pslist)

if plot_evolution is not None:
    # Fill with a value less than any possible value
    zwin = np.full(nfiles, np.finfo(float).min)
    quantity = np.full(nfiles, np.finfo(float).min)

if plot_particle_evolution is not None:
    # Fill with a value less than any possible value
    zbar = np.full(nfiles, np.finfo(float).min)
    xstd = np.full(nfiles, np.finfo(float).min)

# Loop over files, splitting plotfile list among MPI ranks
# - plot field snapshot
# - store window position and field max in arrays
for count, filename in enumerate(file_list):
    if count%size != rank:
        continue

    plot_snapshot( filename )
    if plot_evolution is not None:
        zwin[count], quantity[count] = get_evolution_quantity( filename, plot_evolution )
    if plot_particle_evolution is not None:
        zbar[count], xstd[count] = get_particle_evolution_quantity(filename, plot_particle_evolution)

if plot_evolution is not None:
    zwin, quantity = reduce_evolved_quantity(zwin, quantity)
    if rank == 0:
        plot_evolved_quantity(zwin, quantity)

if plot_particle_evolution is not None:
    zbar, xstd = reduce_evolved_quantity(zbar, xstd)
    if rank == 0:
        plot_particle_evolved_quantity(zbar, xstd)

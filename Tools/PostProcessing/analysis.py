# Main modules
import inspect
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as scc
import re
import yt

# For back-transformed diagnostics
import read_raw_data

from testcase import TestCase

class Analysis:

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------

    def __init__( self, path, backtransf=False, testcase='none' ):
        self.backtransf = backtransf
        if not self.backtransf: # lab frame data from lab frame simulation
               print()
               self.ds = yt.load( path )
               self.header = path + 'Header'
               self.ad = self.ds.covering_grid( level=0, \
                                                left_edge=self.ds.domain_left_edge, \
                                                dims=self.ds.domain_dimensions )
               self.is2D = True if self.ds.domain_dimensions[2] == 1 else False
               self.time = int( (re.search( 'diags/diag(.*)/', path )).group(1) )
               self.extract_domain_2D_edges()
               self.get_fields()
               self.testcase = TestCase( testcase )
        else: # lab frame data from boosted frame simulation
               self.snapshot = path
               self.header = re.sub( (re.search( 'snapshots/(.*)', path )).group(1), 'Header', path )
               self.allrd, self.info = read_raw_data.read_lab_snapshot( self.snapshot, self.header )
               self.is2D = True if self.info['ny'] == 2 else False
               self.time = int( (re.search( 'snapshots/snapshot(.*)/', path )).group(1) )
               self.extract_domain_2D_edges()
               self.get_fields()

    #---------------------------------------------------------------------------
    # Main methods
    #---------------------------------------------------------------------------

    def check_error( self, error, tolerance ):
        print( '  tolerance = %g' %(tolerance) )
        if ( error <= tolerance ):
           print( '\n  relative error <= tolerance: test passed' )
        else:
           print( '\n  relative error > tolerance: test failed' )

    def compute_energy( self ):
        energy_E = 0.5 * np.sum( (self.Ex**2+self.Ey**2+self.Ez**2) ) * scc.epsilon_0
        energy_B = 0.5 * np.sum( (self.Bx**2+self.By**2+self.Bz**2) ) / scc.mu_0
        self.energy = energy_E + energy_B

    def compute_error( self, val, ref ):
        error = np.abs( ref - val ) / np.abs( ref )
        print( "\n  numerical value = %g" %(val) )
        print( "  reference value = %g" %(ref) )
        print( "\n  rel error = %g" %(error) )
        return error

    def compute_reflectivity( self, energy_ref ):
        self.compute_energy()
        self.reflectivity = self.energy / energy_ref

    def extract_domain_2D_edges( self ):
        if not self.backtransf: # lab frame data from lab frame simulation
               if self.is2D: # 2D simulation
                  self.edges = np.array( [ \
                               np.array( [ (self.ds.domain_left_edge [0]).item(), \
                                           (self.ds.domain_right_edge[0]).item(), \
                                           (self.ds.domain_left_edge [1]).item(), \
                                           (self.ds.domain_right_edge[1]).item() ] ) ] )
               else: # 3D simulation
                  self.edges = np.array( [ \
                               np.array( [ (self.ds.domain_left_edge [0]).item(), \
                                           (self.ds.domain_right_edge[0]).item(), \
                                           (self.ds.domain_left_edge [1]).item(), \
                                           (self.ds.domain_right_edge[1]).item() ] ), \
                               np.array( [ (self.ds.domain_left_edge [0]).item(), \
                                           (self.ds.domain_right_edge[0]).item(), \
                                           (self.ds.domain_left_edge [2]).item(), \
                                           (self.ds.domain_right_edge[2]).item() ] ), \
                               np.array( [ (self.ds.domain_left_edge [1]).item(), \
                                           (self.ds.domain_right_edge[1]).item(), \
                                           (self.ds.domain_left_edge [2]).item(), \
                                           (self.ds.domain_right_edge[2]).item() ] ) ] )
        else: # lab frame data from boosted frame simulation
               if self.is2D: # 2D simulation
                  self.edges = np.array( [ \
                               np.array( [ self.info['xmin'], self.info['xmax'], \
                                           self.info['zmin'], self.info['zmax'] ] ) ] )
               else: # 3D simulation
                  self.edges = np.array( [ \
                               np.array( [ self.info['xmin'], self.info['xmax'], \
                                           self.info['ymin'], self.info['ymax'] ] ), \
                               np.array( [ self.info['xmin'], self.info['xmax'], \
                                           self.info['zmin'], self.info['zmax'] ] ), \
                               np.array( [ self.info['ymin'], self.info['ymax'], \
                                           self.info['zmin'], self.info['zmax'] ] ) ] )

    #TODO: loop over variable names
    def get_fields( self ):
        if not self.backtransf: # lab frame data from lab frame simulation
               kl = list( kk[1] for kk in list( self.ds.field_list ) )
               if ( next( ( kk for kk in kl if kk == 'Ex'), None ) ):
                  self.Ex = self.ad['boxlib','Ex'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'Ey'), None ) ):
                  self.Ey = self.ad['boxlib','Ey'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'Ez'), None ) ):
                  self.Ez = self.ad['boxlib','Ez'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'Bx'), None ) ):
                  self.Bx = self.ad['boxlib','Bx'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'By'), None ) ):
                  self.By = self.ad['boxlib','By'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'Bz'), None ) ):
                  self.Bz = self.ad['boxlib','Bz'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'jx'), None ) ):
                  self.Jx = self.ad['boxlib','jx'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'jy'), None ) ):
                  self.Jy = self.ad['boxlib','jy'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'jz'), None ) ):
                  self.Jz = self.ad['boxlib','jz'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'rho'), None ) ):
                  self.rho = self.ad['boxlib','rho'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'divE'), None ) ):
                  self.divE = self.ad['boxlib','divE'].v.squeeze()
               if ( next( ( kk for kk in kl if kk == 'divB'), None ) ):
                  self.divB = self.ad['boxlib','divB'].v.squeeze()
        else: # lab frame data from boosted frame simulation
               kl = list( self.allrd.keys() )
               if ( next( ( kk for kk in kl if kk == 'Ex'), None ) ):
                  self.Ex = self.allrd['Ex']
               if ( next( ( kk for kk in kl if kk == 'Ey'), None ) ):
                  self.Ey = self.allrd['Ey']
               if ( next( ( kk for kk in kl if kk == 'Ez'), None ) ):
                  self.Ez = self.allrd['Ez']
               if ( next( ( kk for kk in kl if kk == 'Bx'), None ) ):
                  self.Bx = self.allrd['Bx']
               if ( next( ( kk for kk in kl if kk == 'By'), None ) ):
                  self.By = self.allrd['By']
               if ( next( ( kk for kk in kl if kk == 'Bz'), None ) ):
                  self.Bz = self.allrd['Bz']
               if ( next( ( kk for kk in kl if kk == 'jx'), None ) ):
                  self.Jx = self.allrd['jx']
               if ( next( ( kk for kk in kl if kk == 'jy'), None ) ):
                  self.Jy = self.allrd['jy']
               if ( next( ( kk for kk in kl if kk == 'jz'), None ) ):
                  self.Jz = self.allrd['jz']
               if ( next( ( kk for kk in kl if kk == 'rho'), None ) ):
                  self.rho = self.allrd['rho']
               if ( next( ( kk for kk in kl if kk == 'divE'), None ) ):
                  self.divE = self.allrd['divE']
               if ( next( ( kk for kk in kl if kk == 'divB'), None ) ):
                  self.divB = self.allrd['divB']

    def perform_analysis( self ):
        if ( self.testcase.name == 'pml_x_psatd' ):
           self.compute_reflectivity( self.testcase.energy_ref )
           error = self.compute_error( self.reflectivity, self.testcase.reflec_ref )
           self.check_error( error, tolerance=5.e-2 )

    #---------------------------------------------------------------------------
    # Plotting
    #---------------------------------------------------------------------------

    def plot_field( self, field, plane='none', label='none', save=False, name='figure.png' ):
        fg = plt.figure()
        ax = fg.add_subplot( 111 )
        vlmin = field.min()
        vlmax = field.max()
        clmin = vlmin #self.compute_clevels_limits( vlmin )
        clmax = vlmax #self.compute_clevels_limits( vlmax )
        clevels = np.linspace( clmin, clmax, 101 )
        if ( plane == 'none' ):
           return
        if ( plane == 'xy' ):
           if self.is2D:
              return
           im = ax.imshow( field.transpose(), extent=self.edges[0], vmin=clmin, vmax=clmax, aspect='auto' )
           ax.set_xlabel( r'$x$' )
           ax.set_ylabel( r'$y$' )
        if ( plane == 'xz' ):
           if self.is2D:
              im = ax.imshow( field.transpose(), extent=self.edges[0], vmin=clmin, vmax=clmax, aspect='auto' )
           else:
              im = ax.imshow( field.transpose(), extent=self.edges[1], vmin=clmin, vmax=clmax, aspect='auto' )
           ax.set_xlabel( r'$x$' )
           ax.set_ylabel( r'$z$' )
        if ( plane == 'yz' ):
           if self.is2D:
              return
           im = ax.imshow( field.transpose(), extent=self.edges[2], vmin=clmin, vmax=clmax, aspect='auto' )
           ax.set_xlabel( r'$y$' )
           ax.set_ylabel( r'$z$' )
        #field_name = retrieve_name( field )[0] # read name string from variable
        if ( label == 'none' ):
           ax.set_title( r'$(t=%d)$: min$=%g$, max$=%g$' %(self.time,vlmin,vlmax), fontsize=10 )
        else:
           ax.set_title( label+r'$(t=%d)$: min$=%g$, max$=%g$' %(self.time,vlmin,vlmax), fontsize=10 )
        ax.xaxis.set_major_locator( plt.MaxNLocator(5) )
        ax.yaxis.set_major_locator( plt.MaxNLocator(5) )
        cb = fg.colorbar( im, ax=ax )
        #cb = fg.colorbar( im, ax=ax, orientation='horizontal', pad=0.2 )
        #tk = np.linspace( clmin, clmax, 10 )
        #cb.set_ticks( tk )
        cb.ax.tick_params()
        cb.ax.yaxis.get_offset_text().set()
        cb.ax.yaxis.offsetText.set_ha( 'center' )
        cb.ax.yaxis.offsetText.set_va( 'bottom' )
        cb.formatter.set_powerlimits( (0,0) )
        cb.update_ticks()
        fg.tight_layout()
        fg.show()
        if save:
           fg.savefig( name, dpi=100 )

    #---------------------------------------------------------------------------
    # Utilities
    #---------------------------------------------------------------------------

    def compute_clevels_limits( self, value ):
        power = np.floor( np.log10( np.abs( value ) ) )
        limit = np.sign( value ) * np.ceil( np.abs( value ) / 10**power ) * 10**power
        return limit

    def retrieve_name( self, var ):
        callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
        return [ var_name for var_name, var_val in callers_local_vars if var_val is var ]

    def set_LaTeX( self ):
        from matplotlib import rc
        rc( 'font', **{'family':'sans-serif','sans-serif':['Helvetica']} )
        rc( 'text', usetex=True )
        mpl.rcParams[ 'text.latex.preamble' ] = [ r"\usepackage{amsmath}" ]

    def set_style( self, fs=10, lw=1.0 ):
        self.fs = fs
        self.lw = lw

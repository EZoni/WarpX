# Main modules
import inspect
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as scc
import re
import sys
import yt

# For back-transformed diagnostics
import read_raw_data

from analysis import Analysis
from analysis import plot_field
from testcase import TestCase

path = []
for a in sys.argv[1:]:
   path.append( a )

backtransf = []
for p in path:
    backtransf.append( True if re.findall( 'lab_frame_data', p ) else False )

data = []
indx = 0
for p in path:
    data.append( Analysis( p, backtransf[indx] ) )
    indx += 1

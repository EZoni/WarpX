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
from testcase import TestCase

path = sys.argv[1]

backtransf = True if re.findall( 'lab_frame_data', path ) else False

this = Analysis( path, backtransf )

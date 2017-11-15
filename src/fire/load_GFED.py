#===============================================================================
# load_GFED.py
#-------------------------------------------------------------------------------
# @author D. T. Milodowski, November 2017
# This is a set of functions to load in GFED4 burned area data into
# numpy arrays. These arrays have three dimensions: lat, long, time.
#
# Data references:
# - Giglio, L., J.T. Randerson, and G.R. van der Werf, (2013), J. Geophys. Res. 
#   Biogeosci., 118, 317328, doi:10.1002/jgrg.20042.
# - Randerson, J.T., Y. Chen, G.R. van derWerf, B.M. Rogers, and D.C. Morton
#   (2012), J. Geophys. Res., 117, G04012, doi:10.1029/2012JG002128.
# - van der Werf, G.R., Randerson, J.T., Giglio, L., et al., (2017), Earth Syst.. 
#   Sci. Data, 9, 697-720, https://doi.org/10.5194/essd-9-697-2017.
#===============================================================================
# import standard libraries
import numpy as np

# import netcdf libraries
from netCDF4 import Dataset


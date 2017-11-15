#===============================================================================
# load_GFED.py
#-------------------------------------------------------------------------------
# @author D. T. Milodowski, November 2017
# This is a set of functions to load in GFED burned area data into
# numpy arrays. These arrays have three dimensions: lat, long, time.
#===============================================================================
# import standard libraries
import numpy as np

# import netcdf libraries
from netCDF4 import Dataset


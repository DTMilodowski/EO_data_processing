#===============================================================================
# load_GFED.py
#-------------------------------------------------------------------------------
# @author D. T. Milodowski, November 2017
# This is a set of functions to load in GFED4 burned area data into numpy arrays 
# These arrays have three dimensions: lat, long, time, with a monthly timestep
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

# Function to load GFED4 monthly data. There are two potential variables here:
# - burned area expressed as fraction of a pixel: BurnedFraction
# - burned area expressed as the area in square metres: BurnedArea
# Files collate annual data, with timestep indicated by the day of the year.
def load_GFED4_monthly(path2files,variable,start_month,start_year,end_month,end_year):

    # first of all obtain the start and end date of the time series
    start_date = np.datetime64('%04i-%02i' % (start_year,start_month))
    end_date = (np.datetime64('%04i-%02i' % (end_year,end_month))+np.timedelta64(1,'M'))


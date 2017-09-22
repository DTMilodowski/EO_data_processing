#===============================================================================
# load_ERAinterim.py
#-------------------------------------------------------------------------------
# @author D. T. Milodowski, September 2017
# This is a set of functions to load in ERAInterim meteorological data into
# numpy arrays. These arrays have three dimensions: lat, long, time.
#===============================================================================
# import standard libraries
import numpy as np

# import netcdf libraries
from netCDF4 import Dataset

#-------------------------------------------------------------------------------
# Script to load in daily ERA interim.
# Variable names permitted are:
# - mx2t
# - mn2t
# - t2m
# - u10w (E-W)
# - v10w (N-s)
# - d2m
# - psurf
# - ssrd
def load_ERAinterim_daily(path2files,variable,start_month,start_year,end_month,end_year):

    month = np.arange(start_month,end_month+1)
    year = np.arange(start_year,end_year+1)
    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,month.size):
            NetCDF_file = path2files+'/'+variable+'_'+str(month[mm]).zfill(2)+str(year[yy])+'.nc'
            dataset = Dataset(NetCDF_file)
            tt+=1

return met_array

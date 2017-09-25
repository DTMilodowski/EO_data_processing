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

    # first of all find the first and last tile, and obtain the start and end
    # date of the time series
    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/%s_%02i%04i.nc' % (path2files, variable,start_month,start_year)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/%s_%02i%04i.nc' % (path2files, variable,end_month,end_year)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg,'h')).astype('datetime64[D]')
    
    # create date list
    year = np.arange(start_year,end_year+1)
    date = np.arange(start_date,end_date)
    # create host array for met data    
    eravar = np.zeros(date.size)*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                NetCDF_file = '%s/%s_%02i%04i.nc' % (path2files, variable,month[mm],year[yy])
                print NetCDF_file
                dataset = Dataset(NetCDF_file)
                tt+=1

    return met_array

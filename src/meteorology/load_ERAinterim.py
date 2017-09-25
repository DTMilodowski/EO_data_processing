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
    metvar = np.zeros(date.size)*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                NetCDF_file = '%s/%s_%02i%04i.nc' % (path2files, variable,month[mm],year[yy])
                print NetCDF_file

                # note that scale and offsets automatically applied when reading in
                # data in this way
                dataset = Dataset(NetCDF_file)
                
                if variable == 'ssrd':
                    N = dataset.variables['time'][:].size/2
                else:
                    N = dataset.variables['time'][:].size/4

                eravar = dataset.variables[variable]
                    
                for ii in range(0,N):
                    # now need to use variable specific processing chain
                    # - minimum and maximum temperatures - four per day, take
                    #   minimum and maximum respectively
                    if variable == 'mn2t':
                        metvar[tt] = np.min(eravar[ii*4:(ii+1)*4])
                    if variable == 'mx2t':
                        metvar[tt] = np.max(eravar[ii*4:(ii+1)*4])
                    # - instantaneous temperatures and dewpoint temperatures,
                    #   wind speed, surface pressure;
                    #   four per day, take average
                    if variable in ['t2m','d2m','u10w','v10w']:
                        metvar[tt] = np.mean(eravar[ii*4:(ii+1)*4])
                    # - ssrd - only two timesteps per day, which need to be summed
                    if variable == 'ssrd':
                        metvar[tt] = np.sum(eravar[ii*2:(ii+1)*2])
                        
                    # iterate timestep
                    tt+=1

    return date,metvar


# Calculate rh based on t2m and d2m
def calculate_rh_daily(path2files,variable,start_month,start_year,end_month,end_year):
    
    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/d2m_%02i%04i.nc' % (path2files,start_month,start_year)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/d2m_%02i%04i.nc' % (path2files,end_month,end_year)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg,'h')).astype('datetime64[D]')
    
    # create date list
    year = np.arange(start_year,end_year+1)
    date = np.arange(start_date,end_date)
    # create host array for met data    
    rh_daily = np.zeros(date.size)*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                d2m_file = '%s/d2m_%02i%04i.nc' % (path2files,month[mm],year[yy])
                t2m_file = '%s/d2m_%02i%04i.nc' % (path2files,month[mm],year[yy])
                # note that scale and offsets automatically applied when reading
                # data in this way
                ds_d2m = Dataset(d2m_file)
                ds_t2m = Dataset(t2m_file)

                N = dataset.variables['time'][:].size/4

                d2m = ds_d2m.variables[variable][:] - 273.15 # convert from K to oC
                t2m = ds_t2m.variables[variable][:] - 273.15 # convert from K to oC
                es_T = 610.94*np.exp((17.625*t2m)/(243.04+t2m))
                es_Td = 610.94*np.exp((17.625*d2m)/(243.04+d2m))
                rh = 100.*(es_Td/es_T)
                for ii in range(0,N):          
                    rh_daily = np.mean(rh[ii*4:(ii+1)*4])
                    tt+=1

    return rh_daily


"""
                    # Also options for variables that need calculating on the fly
                    # - relative humidity
                    if variable == 'rh':
                        es_T = 610.94*exp((17.625*tair_site_unpacked)/(243.04+tair_site_unpacked))
      
                        es_Td = 610.94*exp((17.625*dewpoint_site_unpacked)/(243.04+dewpoint_site_unpacked))
                        
                        rh_site = 100*(es_Td/es_T)
"""

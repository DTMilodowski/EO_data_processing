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
# - prcp
def load_ERAinterim_daily(path2files,variable,start_month,start_year,end_month,end_year):

    # first of all find the first and last tile, and obtain the start and end
    # date of the time series
    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/%s_%04i%02i.nc' % (path2files, variable,start_year,start_month)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/%s_%04i%02i.nc' % (path2files, variable,end_year,end_month)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg,'h')).astype('datetime64[D]')

    # Make exceptions for prcp, ssrd and windspd since these variables are
    # cumulative or average 12 hr values based on previous 12 hours. This means
    # that final timestep is 00:00 which gives an end date a day later than
    # required
    if variable in ['prcp','u10w','v10w','ssrd']:
        end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg-12,'h')).astype('datetime64[D]')

    # create date list
    year = np.arange(start_year,end_year+1)
    date = np.arange(start_date,end_date+np.timedelta64(1,'D'))

    # create host array for met data
    lat = np.asarray(dataset.variables['latitude'])
    lon = np.asarray(dataset.variables['longitude'])
    lon[lon>180]-=360

    metvar = np.zeros((date.size,lat.size,lon.size))*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                NetCDF_file = '%s/%s_%04i%02i.nc' % (path2files, variable,year[yy],mm+1)
                #print(NetCDF_file)

                # note that scale and offsets automatically applied when reading in
                # data in this way
                dataset = Dataset(NetCDF_file)

                if variable in ['ssrd','u10w','v10w','mn2t','mx2t','prcp']:
                    N = dataset.variables['time'][:].size/2
                else:
                    N = dataset.variables['time'][:].size/4

                varcode = variable
                if variable == 'prcp':
                    varcode = 'tp'
                elif variable == 'u10w':
                    varcode = 'u10'
                elif variable == 'v10w':
                    varcode = 'v10'
                elif variable == 'psurf':
                    varcode = 'sp'
                eravar = dataset.variables[varcode]

                for ii in range(0,N):
                    # now need to use variable specific processing chain
                    # - minimum and maximum temperatures - two per day, take
                    #   minimum and maximum respectively
                    if variable == 'mn2t':
                        metvar[tt,:,:] = np.min(eravar[ii*2:(ii+1)*2,:,:],axis=0)
                    if variable == 'mx2t':
                        metvar[tt,:,:] = np.max(eravar[ii*2:(ii+1)*2,:,:],axis=0)
                    # - instantaneous temperatures and dewpoint temperatures,
                    #   surface pressure;
                    #   four per day, take average
                    if variable in ['t2m','d2m','psurf']:
                        metvar[tt,:,:] = np.mean(eravar[ii*4:(ii+1)*4,:,:],axis=0)
                    # - ssrd - only two timesteps per day, which need to be summed
                    if variable in ['prcp','ssrd']:
                        metvar[tt,:,:] = np.sum(eravar[ii*2:(ii+1)*2,:,:],axis=0)
                    # - wind speeds - only two timesteps per day, average
                    if variable in ['u10w','v10w']:
                        metvar[tt,:,:] = np.mean(eravar[ii*2:(ii+1)*2,:,:],axis=0)

                    # iterate timestep
                    tt+=1

    if variable in ['t2m','d2m','mn2t','mx2t']:
        metvar -= 273.15 # convert from K to oC
    return date,lat,lon,metvar


# Calculate rh based on t2m and d2m
def calculate_rh_daily(path2files,start_month,start_year,end_month,end_year):

    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/d2m_%04i%02i.nc' % (path2files,start_year,start_month)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/d2m_%04i%02i.nc' % (path2files,end_year,end_month)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg,'h')).astype('datetime64[D]')

    # create date list
    year = np.arange(start_year,end_year+1)
    date = np.arange(start_date,end_date+np.timedelta64(1,'D'))
    # create host array for met data
    lat = np.asarray(dataset.variables['latitude'])
    lon = np.asarray(dataset.variables['longitude'])
    lon[lon>180]-=360

    rh_daily = np.zeros((date.size,lat.size,lon.size))*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                d2m_file = '%s/d2m_%04i%02i.nc' % (path2files,year[yy],mm+1)
                t2m_file = '%s/t2m_%04i%02i.nc' % (path2files,year[yy],mm+1)
                # note that scale and offsets automatically applied when reading
                # data in this way
                ds_d2m = Dataset(d2m_file)
                ds_t2m = Dataset(t2m_file)

                N = ds_d2m.variables['time'][:].size/4

                d2m = ds_d2m.variables['d2m'][:] - 273.15 # convert from K to oC
                t2m = ds_t2m.variables['t2m'][:] - 273.15 # convert from K to oC
                es_T = 610.94*np.exp((17.625*t2m)/(243.04+t2m))
                es_Td = 610.94*np.exp((17.625*d2m)/(243.04+d2m))
                rh = 100.*(es_Td/es_T)
                for ii in range(0,N):
                    rh_daily[tt,:,:] = np.mean(rh[ii*4:(ii+1)*4,:,:],axis=0)
                    tt+=1

    rh_daily[rh_daily>100]=100.
    return date,lat,lon, rh_daily


# Calculate vpd based on t2m and d2m
def calculate_vpd_daily(path2files,start_month,start_year,end_month,end_year):

    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/d2m_%04i%02i.nc' % (path2files,start_year,start_month)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/d2m_%04i%02i.nc' % (path2files,end_year,end_month)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg,'h')).astype('datetime64[D]')

    # create date list
    year = np.arange(start_year,end_year+1)
    date = np.arange(start_date,end_date+np.timedelta64(1,'D'))
    # create host array for met data
    lat = np.asarray(dataset.variables['latitude'])
    lon = np.asarray(dataset.variables['longitude'])
    lon[lon>180]-=360

    vpd_daily = np.zeros((date.size,lat.size,lon.size))*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                d2m_file = '%s/d2m_%04i%02i.nc' % (path2files,year[yy],mm+1)
                t2m_file = '%s/t2m_%04i%02i.nc' % (path2files,year[yy],mm+1)
                # note that scale and offsets automatically applied when reading
                # data in this way
                ds_d2m = Dataset(d2m_file)
                ds_t2m = Dataset(t2m_file)

                N = int(ds_d2m.variables['time'][:].size/4)

                d2m = ds_d2m.variables['d2m'][:] - 273.15 # convert from K to oC
                t2m = ds_t2m.variables['t2m'][:] - 273.15 # convert from K to oC

                es = 610.94*np.exp((17.625*t2m)/(243.04+t2m))
                ea = 610.94*np.exp((17.625*d2m)/(243.04+d2m))
                vpd = es-ea
                for ii in range(0,N):
                    vpd_daily[tt,:,:] = np.mean(vpd[ii*4:(ii+1)*4,:,:],axis=0)
                    tt+=1

    return date,lat,lon, vpd_daily

# Calculate vpd based on t2m and d2m
def calculate_wind_speed_daily(path2files,start_month,start_year,end_month,end_year):

    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/u10w_%04i%02i.nc' % (path2files,start_year,start_month)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/u10w_%04i%02i.nc' % (path2files,end_year,end_month)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg-12,'h')).astype('datetime64[D]')

    # create date list
    year = np.arange(start_year,end_year+1)
    date = np.arange(start_date,end_date+np.timedelta64(1,'D'))
    # create host array for met data
    lat = np.asarray(dataset.variables['latitude'])
    lon = np.asarray(dataset.variables['longitude'])
    lon[lon>180]-=360

    w_daily = np.zeros((date.size,lat.size,lon.size))*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                u_file = '%s/u10w_%04i%02i.nc' % (path2files,year[yy],mm+1)
                v_file = '%s/v10w_%04i%02i.nc' % (path2files,year[yy],mm+1)
                # note that scale and offsets automatically applied when reading
                # data in this way
                ds_u = Dataset(u_file)
                ds_v = Dataset(v_file)

                N = int(ds_u.variables['time'][:].size/2)

                u = ds_u.variables['u10'][:]
                v = ds_v.variables['v10'][:]

                w = np.sqrt(u*u + v*v)
                for ii in range(0,N):
                    w_daily[tt,:,:] = np.mean(w[ii*2:(ii+1)*2,:,:],axis=0)
                    tt+=1

    return date,lat,lon, w_daily

def load_trmm_3B42_daily(path2files,start_year,end_year):

    start_date = np.datetime64('%04i-01-01' % start_year)
    end_date = np.datetime64('%04i-01-01' % (end_year+1))

    # create date list
    year = np.arange(start_year,end_year+1)
    date = np.arange(start_date,end_date)
    # create host array for met data

    pptn_file = '%s3B42_Daily.%04i%02i%02i.7.nc4' % (path2files,start_year,1,1)
    ds = Dataset(pptn_file)
    lat = np.asarray(ds.variables['lat'])
    lon = np.asarray(ds.variables['lon'])
    lon[lon>180]-=360

    pptn = np.zeros((date.size,lat.size,lon.size))*np.nan
    tt = 0
    for yy in range(start_year,end_year+1):
        for mm in range(1,13):
            days = 30
            if mm in [1,3,5,7,8,10,12]:
                days = 31
            elif mm == 2:
                if yy%4 == 0:
                    days = 29
                else:
                    days = 28
            for dd in range(1,days+1):
                pptn_file = '%s3B42_Daily.%04i%02i%02i.7.nc4' % (path2files,yy,mm,dd)
                ds = Dataset(pptn_file)
                pptn[tt] = np.transpose(ds.variables['precipitation'][:])
                tt+=1


    return date,lat,lon, pptn



# Load ERA Interim data for nearest location
def load_nearest_point_ERAInterim_daily(variable,lat,lon,start_date_str,end_date_str,path2files=''):
    if len(path2files)==0:
        # default file path
        path2files='/disk/scratch/local.2/dmilodow/ERAinterim/source_files/0.25deg_Indonesia/'

    start_year,start_month,start_day = start_date_str.split('-')
    end_year,end_month,end_day = end_date_str.split('-')

    # first of all find the first and last tile, and obtain the start and end
    # date of the time series
    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/%s_%s%s.nc' % (path2files, variable,start_year,start_month)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/%s_%s%s.nc' % (path2files, variable,end_year,end_month)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg,'h')).astype('datetime64[D]')

    # Make exceptions for prcp, ssrd and windspd since these variables are
    # cumulative or average 12 hr values based on previous 12 hours. This means
    # that final timestep is 00:00 which gives an end date a day later than
    # required
    if variable in ['prcp','u10w','v10w','ssrd']:
        end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg-12,'h')).astype('datetime64[D]')

    # create date list
    year = np.arange(int(start_year),int(end_year)+1)
    date = np.arange(start_date,end_date+np.timedelta64(1,'D'))

    # create host array for met data
    lats = np.asarray(dataset.variables['latitude'])
    lons = np.asarray(dataset.variables['longitude'])
    lons[lons>180]-=360
    lat_idx = np.abs(lats - lat).argmin()
    lon_idx = np.abs(lons - lon).argmin()

    metvar = np.zeros((date.size))*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt < date.size:
                NetCDF_file = '%s/%s_%04i%02i.nc' % (path2files, variable,year[yy],mm+1)
                #print(NetCDF_file)

                # note that scale and offsets automatically applied when reading in
                # data in this way
                dataset = Dataset(NetCDF_file)

                if variable in ['ssrd','u10w','v10w','mn2t','mx2t','prcp']:
                    N = dataset.variables['time'][:].size/2
                else:
                    N = dataset.variables['time'][:].size/4

                varcode = variable
                if variable == 'prcp':
                    varcode = 'tp'
                elif variable == 'u10w':
                    varcode = 'u10'
                elif variable == 'v10w':
                    varcode = 'v10'
                elif variable == 'psurf':
                    varcode = 'sp'
                eravar = dataset.variables[varcode]

                for ii in range(0,N):
                    # now need to use variable specific processing chain
                    # - minimum and maximum temperatures - two per day, take
                    #   minimum and maximum respectively
                    if variable == 'mn2t':
                        metvar[tt] = np.min(eravar[ii*2:(ii+1)*2,lat_idx,lon_idx])
                    if variable == 'mx2t':
                        metvar[tt] = np.max(eravar[ii*2:(ii+1)*2,lat_idx,lon_idx])
                    # - instantaneous temperatures and dewpoint temperatures,
                    #   surface pressure;
                    #   four per day, take average
                    if variable in ['t2m','d2m','psurf']:
                        metvar[tt] = np.mean(eravar[ii*4:(ii+1)*4,lat_idx,lon_idx])
                    # - ssrd - only two timesteps per day, which need to be summed
                    if variable in ['prcp','ssrd']:
                        metvar[tt] = np.sum(eravar[ii*2:(ii+1)*2,lat_idx,lon_idx])
                    # - wind speeds - only two timesteps per day, average
                    if variable in ['u10w','v10w']:
                        metvar[tt] = np.mean(eravar[ii*2:(ii+1)*2,lat_idx,lon_idx])

                    # iterate timestep
                    tt+=1

    if variable in ['t2m','d2m','mn2t','mx2t']:
        metvar -= 273.15 # convert from K to oC
    return date,metvar

# as above but without the daily aggregation
def load_nearest_point_ERAInterim_raw(variable,lat,lon,start_date_str,end_date_str,path2files=''):
    if len(path2files)==0:
        # default file path
        path2files='/disk/scratch/local.2/dmilodow/ERAinterim/source_files/0.25deg_Indonesia/'

    start_year,start_month,start_day = start_date_str.split('-')
    end_year,end_month,end_day = end_date_str.split('-')

    # first of all find the first and last tile, and obtain the start and end
    # date of the time series
    # Note that ERA Interim times are in gregorian calendar format, hours after
    # 1900-01-01 00:00
    NetCDF_file = '%s/%s_%s%s.nc' % (path2files, variable,start_year,start_month)
    dataset = Dataset(NetCDF_file)
    start_greg = int(dataset.variables['time'][0])
    start_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(start_greg,'h')).astype('datetime64[D]')

    NetCDF_file = '%s/%s_%s%s.nc' % (path2files, variable,end_year,end_month)
    dataset = Dataset(NetCDF_file)
    end_greg = int(dataset.variables['time'][-1])
    end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg,'h')).astype('datetime64[D]')

    # Make exceptions for prcp, ssrd and windspd since these variables are
    # cumulative or average 12 hr values based on previous 12 hours. This means
    # that final timestep is 00:00 which gives an end date a day later than
    # required
    if variable in ['prcp','u10w','v10w','ssrd']:
        end_date = (np.datetime64('1900-01-01 00:00') + np.timedelta64(end_greg-12,'h')).astype('datetime64[D]')

    # create date list
    year = np.arange(int(start_year),int(end_year)+1)
    date = np.arange(start_date,end_date+np.timedelta64(1,'D'))

    # create host array for met data
    lats = np.asarray(dataset.variables['latitude'])
    lons = np.asarray(dataset.variables['longitude'])
    lons[lons>180]-=360
    lat_idx = np.abs(lats - lat).argmin()
    lon_idx = np.abs(lons - lon).argmin()

    # setup array to host the metdata
    if variable in ['ssrd','u10w','v10w','mn2t','mx2t','prcp']:
        Ntot=date.size*2
    else:
        Ntot=date.size*4
    metvar = np.zeros((Ntot))*np.nan

    tt = 0
    for yy in range(0,year.size):
        for mm in range(0,12):
            if tt >= Ntot:
                break
            else:
                NetCDF_file = '%s/%s_%04i%02i.nc' % (path2files, variable,year[yy],mm+1)
                #print(NetCDF_file)

                # note that scale and offsets automatically applied when reading in
                # data in this way
                dataset = Dataset(NetCDF_file)
                N = dataset.variables['time'][:].size

                varcode = variable
                if variable == 'prcp':
                    varcode = 'tp'
                elif variable == 'u10w':
                    varcode = 'u10'
                elif variable == 'v10w':
                    varcode = 'v10'
                elif variable == 'psurf':
                    varcode = 'sp'
                eravar = dataset.variables[varcode]
                metvar[tt:tt+N]=eravar[:,lat_idx,lon_idx]
                # iterate timestep
                tt+=N

    if variable in ['t2m','d2m','mn2t','mx2t']:
        metvar -= 273.15 # convert from K to oC
    return metvar
"""
# A file to clip ERA Interim to a specified bounding box
def clip_netcdf_to_bbox(netcdf_in,netcdf_out,N,S,E,W,lat_var = 'latitude',lon_var = 'longitude',time_var = 'time'):

    src = Dataset(netcdf_in)
    dst = Dataset(netcdf_out, "w")

    # copy attributes
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    # copy dimensions
    for name, dimension in src.dimensions.iteritems():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited else None))

    # copy all file data except for the excluded
    for name, variable in src.variables.iteritems():
        print name, variable
        if name =
        if name not in toexclude:
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst.variables[name][:] = src.variables[name][:]
"""

# This function regrids daily TRMM precipitation data to a specified spatial and temporal resolution

# import required python libraries
import numpy as np
import sys
import time
from netCDF4 import Dataset
# import my own libraries
import data_io as io
import auxilliary_functions as aux
import geospatial_utility_tools as geo

# Where is the data?
datadir = '/disk/scratch/local.2/TRMM/source_files/'

# where will we store data?
savedir = '/disk/scratch/local.2/dmilodow/TRMM/regridded/'

# what will the id tag be for the file prefix
prefix = 'TRMM_weekly_Mexico_0.125deg'

# What is the bounding box of ROI? In this case use Mexico
N = 33.
S = 14.
E = -86.
W = -118.

# Time series info
start = np.datetime64('2001-01-01','D')
end = np.datetime64('2016-01-01','D')
step = np.timedelta64(7,'D')
steps  = np.arange(start,end,step)
n_steps = steps.size
weeks = np.arange(n_steps)
n_days = 7

# Create the host array
dY = 0.125
dX = 0.125
lat_host = np.arange(S,N,dY)+dY/2. # shifting to cell centre
lon_host = np.arange(W,E,dX)+dX/2. # shifting to cell centre

rows_host = lat_host.size
cols_host = lon_host.size
TRMM_regrid = np.zeros((rows_host,cols_host,n_steps))

# loop through tsteps
for ss in range(0, n_steps):
    # clip time series for last step
    if ss+1 == n_steps:
        n_days=(end-steps[-1]).astype('int')

    # load the datafile
    for dd in range(0,n_days):
        date = steps[ss]+np.timedelta64(dd,'D')
        year = date.astype('datetime64[Y]').astype(int) + 1970
        month = date.astype('datetime64[M]').astype(int) % 12 + 1
        day = (date - date.astype('datetime64[M]')).astype(int) + 1
        tile = datadir+'3B42_Daily.'+str(year).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+'.7.nc4'
        print tile
        pptn, lat, lon  = io.load_TRMM_NetCDF(tile)

        # clip TRMM tile to extent and interpolate to finer resolution using bilinear interpolation
        x_pts = (lon_host - np.min(lon))/np.abs(lon[1]-lon[0])
        y_pts = (lat_host - np.min(lat))/np.abs(lat[1]-lat[0])

        x,y = np.meshgrid(x_pts,y_pts)

        TRMM_regrid[:,:,ss] += geo.bilinear_interpolate(pptn,x,y)


np.savez(TRMM_regrid,savedir+'TRMM_Mexico_regrid_2001_2015.npz')

#---------------------------------------------------------------------------------------------------------------------------
# Write to netcdf file
if ('%s.nc' % (savedir+prefix)) in os.listdir(os.getcwd()):
    os.remove('%s.nc' % (savedir+prefix))

fnc=Dataset('%s.nc' % savedir+prefix,'w')

fnc.createDimension('latitude',lat_host.shape[0])
fnc.createDimension('longitude',lon_host.shape[0])
fnc.createDimension('time',len(weeks))

fnc.createVariable('latitude','d',dimensions=['latitude'])
fnc.variables['latitude'][:]=lat_host
fnc.variables['latitude'].long_name='Latitude N'
fnc.variables['latitude'].units='degrees'

fnc.createVariable('longitude','d',dimensions=['longitude'])
fnc.variables['longitude'][:]=lon_host
fnc.variables['longitude'].long_name='Longitude E'
fnc.variables['longitude'].units='degrees'

fnc.createVariable('time','i',dimensions=['time'])
fnc.variables['time'][:]=weeks
fnc.variables['time'].units='weeks'
fnc.variables['time'].long_name='number of weeks that have elapsed since 2001-01-01'

fnc.createVariable('pptn','d',dimensions=['time','latitude','longitude'], zlib = True, complevel = 1)
fnc.variables['pptn'][:,:,:]=TRMM_regrid
fnc.variables['pptn'].long_name='Total weekly precipitation'
fnc.variables['pptn'].units='fraction of grid cell in which forest loss has occurred within that month'
fnc.variables['pptn'].missing_value=0.

fnc.production_date='%s %s' % (time.asctime(),time.tzname[0])
fnc.production_software='Python %s - netCDF4 library %s' % (sys.version,netCDF4.getlibversion())
fnc.production_source='TRMM'
fnc.production_method='This file contains the total weekly precipitation derived from daily TRMM estimates. Data are regridded to a %4.2fx%4.2f grid. % (dX,dY)

fnc.sync()
fnc.close()
print "DONE!"

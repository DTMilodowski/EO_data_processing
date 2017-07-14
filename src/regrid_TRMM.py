# This function regrids daily TRMM precipitation data to a specified spatial and temporal resolution

# import required python libraries 
import numpy as np
import sys

# import my own libraries
import data_io as io
import auxilliary_functions as aux
import geospatial_utility_tools as geo

# Where is the data?
datadir = '/disk/scratch/local.2/dmilodow/TRMM/source_files/'

# where will we store data?
savedir = '/disk/scratch/local.2/dmilodow/TRMM/regridded/'

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
n_days = 7

# Create the host array
dY = 0.125
dX = 0.125
lat_host = np.arange(S,N,dY)+dY/2. # shifting to cell centre
lon_host = np.arange(W,E,dX)+dX/2. # shifting to cell centre
#areas_host = geo.calculate_cell_area_array(lat_host,lon_host, area_scalar = 1./10.**6,cell_centred=True)

rows_host = lat_host.size
cols_host = lon_host.size
TRMM_regrid = np.zeros((rows_host,cols_host,n_steps))

# loop through tsteps
for ss in range(0, n_steps):

    # clip time series for last step
    if ss+1 == n_steps:
        n_days=int(end-steps[-1])

    # load the datafile
    for dd in range(0,n_days):
        date = steps[ss]+np.timedelta64(dd,'D')

        year = date.astype('datetime64[Y]').astype(int) + 1970
        month = date.astype('datetime64[M]').astype(int) % 12 + 1
        day = (date - date.astype('datetime64[M]')).astype(int) + 1
        
        tile = datadir+'3B42_Daily.'+str(year).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+'.7.nc4'
        pptn, lat, lon  = io.load_TRMM_NetCDF(tile)

        # clip TRMM tile to extent and interpolate to finer resolution using bilinear interpolation
        x_pts = (lon_host - np.min(lon))/np.abs(lon[1]-lon[0])
        y_pts = (lat_host - np.min(lat))/np.abs(lat[1]-lat[0])

        x,y = np.meshgrid(x_pts,y_pts)

        TRMM_regrid[:,:,ss] += geo.bilinear_interpolate(pptn,x,y)

        plt.imshow(TRMM_regrid[:,:,ss],origin='lower');plt.show()

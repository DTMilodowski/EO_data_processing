# This function analyses the global forest watch forest loss dataset to provide
# forest loss time series at a given spatial resolution.  Forest loss is
# calculated based on the fractional area of each pixel for which forest loss
# occurs in a specified timestep.  Initial estimates will be annual, but
# subsequently will move onto monthly estimates using FORMA to distribute the
# change through the year.

# import required python libraries 
import numpy as np
import sys,os
from scipy.ndimage.filters import uniform_filter

import time
from netCDF4 import Dataset
import netCDF4

# import my own libraries
import data_io as io
import auxilliary_functions as aux
import geospatial_utility_tools as geo

# Where is the data?
datadir = '/disk/scratch/local.2/dmilodow/GFW/'

# where will we store data?
savedir = '/disk/scratch/local.2/dmilodow/GFW/regridded/'

# what will the id tag be for the file prefix
prefix = 'GFW_monthly_Mexico_0.125deg'

# What is the bounding box of ROI? In this case use Mexico
N = 33.
S = 14.
E = -86.
W = -118.

# What is the tree cover thrshold for forest loss?
min_treecover = 30 # %

# Time series info
start = 2001
end = 2015
years = np.arange(start,end)
n_years = 2015-2001

# Create the host array
dY = 0.125
dX = 0.125
lat_host = np.arange(S,N,dY)+dY/2. # shifting to cell centre
long_host = np.arange(W,E,dX)+dX/2. # shifting to cell centre
areas_host = geo.calculate_cell_area_array(lat_host,long_host, area_scalar = 1./10.**6,cell_centred=True)

rows_host = lat_host.size
cols_host = long_host.size
regrid = np.zeros((n_years,rows_host,cols_host))

# Read in the list of available GFW tiles
tile_list=np.genfromtxt(datadir+'tile_list.txt',dtype='string')
n_tiles = len(tile_list)

# loop through the tiles
for tt in range(0,n_tiles):
    # only analyse tile of it falls within the ROI
    if io.is_GeoTIFF_within_bbox(datadir+tile_list[tt],N,S,W,E):
        print tile_list[tt]

        # load forest loss year - we read this directly from the tile list
        lossyear_i, geoTrans_i, coord_sys = io.load_GeoTIFF_band_and_georeferencing(datadir+tile_list[tt],band_number=1)

        # clip GFW tile to extent
        print "\t-clipping to bbox extent..."
        lossyear, geoTrans = geo.clip_array_to_bbox(lossyear_i,geoTrans_i,N,S,W,E)

        lossyear_i=None
        geoTrans_i=None

        # calculate cell areas for geographic coordinate system
        print "\t-calculating cell areas..."
        rows,cols=lossyear.shape
        latitude = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3]-0.000001*geoTrans[5],geoTrans[5])+geoTrans[5]/2. # shifting to cell centre
        longitude = np.arange(geoTrans[0],cols*geoTrans[1]+geoTrans[0]-0.000001*geoTrans[1],geoTrans[1])+geoTrans[1]/2. # shifting to cell centre
        # note in the above I added an arbitrarily small alteration to the upper limits to account for rare propagation of float rounding errors
        # This becomes an issue with big datasets like GFW (40000 x 40000 pixels per tile)
 
        areas = np.empty((latitude.size,1))
        areas[:,0] = geo.calculate_cell_area_column(latitude,geoTrans[1], area_scalar = 1./10.**6,cell_centred=True)

        print "\t-finding nearest neighbour..."
        #assign closest point in regrid lat to orig
        closest_lat=np.zeros(rows).astype("int")
        closest_long=np.zeros(cols).astype("int")
        
        for ii,val in enumerate(latitude):
            closest_lat[ii]=np.argsort(np.abs(val-lat_host))[0]
        for jj,val in enumerate(longitude):
            closest_long[jj]=np.argsort(np.abs(val-long_host))[0]

        print "\t-regridding..."
        """
        for ii in range(0,rows_host):
            lat_mask = closest_lat==ii
            if lat_mask.sum()>0:
                for jj in range(0,cols_host):
                    long_mask = closest_long==jj
                    if long_mask.sum() > 0:
                        regrid[ii,jj] += np.sum(lossyear[np.ix_(lat_mask,long_mask)])
        """
        # loop through years and assign change to year - expressed as fraction of pixel deforested in regridded dataset
        for yy in range(0,n_years):
            lossarea = np.zeros((rows,cols))
            lossarea[lossyear==(years[yy]-2000)] = 1.
            lossarea=np.multiply(lossarea,areas)

            for ii in range(0,rows_host):
                lat_mask = closest_lat==ii
                if lat_mask.sum()>0:
                    for jj in range(0,cols_host):
                        long_mask = closest_long==jj
                        if long_mask.sum() > 0:
                            regrid[yy,ii,jj] += np.sum(lossarea[np.ix_(lat_mask,long_mask)])
        
# normalise forest loss to give fraction loss
for yy in range(0,n_years):
    regrid[yy,:,:]/=areas_host

np.savez('GFW_annual.npz',regrid)



regrid = np.load('GFW_annual.npz')['arr_0']
#----------------------------------------------------------------------------------------------------------------------------------------------
# Now load in FORMA.
FORMAfile = '/home/dmilodow/DataStore_GCEL/FORMA/forma-1.0-2005-12-19-2015-08-13.csv'
months,lat_FORMA,lon_FORMA, monthly_degrad = io.grid_FORMA_monthly(FORMAfile,dX,N,S,E,W,start_date = '2006-01-01', end_date='2015-01-01')

FORMA_seasonal = np.zeros((12,monthly_degrad.shape[1],monthly_degrad.shape[2]))
month = (months-months.astype('datetime64[Y]')).astype('int')
# smooth seasonal signal using a moving window that is approx 1 degree by 1 degree
filter_window = np.int(2./dY)
total_FORMA=np.zeros(12)
if filter_window % 2. != 0:
    filter_window+=1
for mm in range(0,12):
    FORMA_seasonal[mm,:,:] = uniform_filter(np.mean(monthly_degrad[month==mm,:,:],axis=0), size=filter_window)
    total_FORMA[mm] = np.mean(monthly_degrad[month==mm,:,:],axis=0).sum()

norm_FORMA=total_FORMA/total_FORMA.sum()
# normalise FORMA
FORMA_sum = np.sum(FORMA_seasonal,axis=0)
FORMA_norm = FORMA_sum.copy()
FORMA_norm[FORMA_sum==0]=1. # note that this just accounts for cases where there is no detected deforestation
FORMA_seasonal_norm = np.zeros(FORMA_seasonal.shape)
for mm in range(0,12):
    FORMA_seasonal_norm[mm,:,:] = np.flipud(FORMA_seasonal[mm,:,:]/FORMA_norm)
    #FORMA_seasonal_norm[mm,:,:] = np.flipud(FORMA_seasonal[mm,:,:])
#----------------------------------------------------------------------------------------------------------------------------------------------
# Now downsample GFW to monthly based on FORMA signal
n_years = regrid.shape[0]
GFW_monthly = np.zeros((n_years*12,regrid.shape[1],regrid.shape[2])) 
month = 0
year = 0
for mm in range(0,n_years*12):
    #GFW_monthly[mm,:,:]=regrid[year,:,:]*FORMA_seasonal_norm[month,:,:]
    GFW_monthly[mm,:,:]=regrid[year,:,:]*norm_FORMA[month]
    month+=1
    if month==12:
        year+=1
        month=0

months_out = np.arange(GFW_monthly.shape[0])
np.savez('GFW_monthly.npz',GFW_monthly)
#----------------------------------------------------------------------------------------------------------------------------------------------
# Now write to netcdf file
if ('%s' % prefix+'.nc') in os.listdir(os.getcwd()):
    os.remove('%s' % savedir+prefix+'.nc')

fnc=Dataset('%s' % savedir+prefix+'.nc','w')

fnc.createDimension('latitude',lat_host.shape[0])
fnc.createDimension('longitude',long_host.shape[0])
fnc.createDimension('time',len(months_out))

fnc.createVariable('latitude','d',dimensions=['latitude'])
fnc.variables['latitude'][:]=lat_host
fnc.variables['latitude'].long_name='Latitude N'
fnc.variables['latitude'].units='degrees'

fnc.createVariable('longitude','d',dimensions=['longitude'])
fnc.variables['longitude'][:]=long_host
fnc.variables['longitude'].long_name='Longitude E'
fnc.variables['longitude'].units='degrees'
    
fnc.createVariable('time','i',dimensions=['time'])
fnc.variables['time'][:]=months_out
fnc.variables['time'].units='months'
fnc.variables['time'].long_name='number of months that have elapsed since 2001-01-01'    

fnc.createVariable('ForestLoss','d',dimensions=['time','latitude','longitude'], zlib = True, complevel = 1)
fnc.variables['ForestLoss'][:,:,:]=GFW_monthly
fnc.variables['ForestLoss'].long_name='GFL/FORMA fusion monthly forest loss'
fnc.variables['ForestLoss'].units='fraction of grid cell in which forest loss has occurred within that month'
fnc.variables['ForestLoss'].missing_value=0.

fnc.production_date='%s %s' % (time.asctime(),time.tzname[0])
fnc.production_software='Python %s - netCDF4 library %s' % (sys.version,netCDF4.getlibversion())
fnc.production_source='UMD Global Forest Loss Dataset & FORMA'
fnc.production_method='This file contains a forest loss product created by fusing Hansen and FORMA.   Data are regridded to a %4.2fx%4.2f grid.  Monthly forest loss are estimated by scaling annual GFL forest loss extents by ratio of FORMA disturbance for each month divided by annual FORMA disturbance in that year' % (dX,dY)

fnc.sync()
fnc.close()

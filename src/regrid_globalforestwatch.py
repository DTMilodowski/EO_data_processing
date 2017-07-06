# This function analyses the global forest watch forest loss dataset to provide
# forest loss time series at a given spatial resolution.  Forest loss is
# calculated based on the fractional area of each pixel for which forest loss
# occurs in a specified timestep.  Initial estimates will be annual, but
# subsequently will move onto monthly estimates using FORMA to distribute the
# change through the year.

# import required python libraries 
import numpy as np
import sys

# import my own libraries
import data_io as io
import auxilliary_functions as aux
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomass/src')
import geospatial_utility_tools as geo

# Where is the data?
datadir = '/disk/scratch/local.2/dmilodow/GFW/'

# where will we store data?
savedir = '/disk/scratch/local.2/dmilodow/GFW/regridded/'

# What is the bounding box of ROI? In this case use Mexico
N = 33.
S = 14.
E = -118.
W = -86.

# What is the tree cover thrshold for forest loss?
min_treecover = 30 # %

# Time series info
start = 2001
end = 2016
years = np.arange(start,end)
n_years = 2016-2001

# Create the host array
dY = 0.175
dX = 0.175
rows_host = int((N-S)/dY)
cols_host = int((W-E)/dX)
regrid = np.zeros((rows_host,cols_host,n_years))
lat_host = np.linspace(S,N,dY)+dY/2. # shifting to cell centre
long_host = np.linspace(W,E,dX)+dX/2. # shifting to cell centre
areas_host = geo.calculate_cell_area_array(lat_host,long_host, area_scalar = 1./10.**6,cell_centred=True)

# Read in the list of available GFW tiles
tile_list=np.genfromtxt(datadir+'tile_list.txt',dtype='string')
n_tiles = len(tile_list)
for ii in range(0,n_tiles):
    # only analyse tile of it falls within the ROI
    if io.is_GeoTIFF_within_bbox(datadir+tile_list[ii],N,S,W,E):
        print tile_list[ii]

        # load forest loss year - we read this directly from the tile list
        lossyear, geoTrans, coord_sys = io.load_GeoTIFF_band_and_georeferencing(datadir+tile_list[ii],band_number=1)
        
        # calculate cell areas for geographic coordinate system
        rows,cols=lossyear.shape
        latitude = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3],geoTrans[5])+geoTrans[5]/2. # shifting to cell centre
        longitude =  np.arange(geoTrans[0],cols*geoTrans[1]+geoTrans[0],geoTrans[1])+geoTrans[1]/2. # shifting to cell centre
        areas = geo.calculate_cell_area_array(latitude,longitude, area_scalar = 1./10.**6,cell_centred=True)

        #assign closest point in regrid lat to orig
        closest_lat=np.zeros(rows).astype("int")
        closest_long=np.zeros(cols).astype("int")
        
        for ii,val in enumerate(lat_host):
            closest_lat[ii]=np.argsort(np.abs(val-regridlat))[0]
        for jj,val in enumerate(long_host):
            closest_lon[jj]=np.argsort(np.abs(val-regridlon))[0]

        # loop through years and assign change to year - expressed as fraction of pixel deforested in regridded dataset
        for yy in range(0,n_years):
            lossarea = areas.copy()
            lossarea[lossyear!=(years[yy]-2000)]=0.
            for ii,lat_ii in enumerate(lat_host):
                lat_mask = closest_lat==lat_ii
                for jj,long_jj in enumerate(long_host):
                    long_mask = closest_long==long_jj
                    regrid[ii,jj,yy] = np.sum(lossarea[lat_mask,long_mask])/areas_host[ii,jj]
            
        

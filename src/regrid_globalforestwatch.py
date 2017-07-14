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
import geospatial_utility_tools as geo

# Where is the data?
datadir = '/disk/scratch/local.2/dmilodow/GFW/'

# where will we store data?
savedir = '/disk/scratch/local.2/dmilodow/GFW/regridded/'

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
regrid = np.zeros((rows_host,cols_host,n_years))

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
                            regrid[ii,jj,yy] += np.sum(lossarea[np.ix_(lat_mask,long_mask)])
        
# normalise forest loss to give fraction loss
for yy in range(0,n_years):
    regrid[:,:,yy]/=areas_host

np.savez('regridded_data',regrid)

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

# Where is the data?
datadir = '/disk/scratch/local.2/dmilodow/GFW/'

# where will we store data?
savedir = '/disk/scratch/local.2/dmilodow/GFW/regridded/'

# What is the bounding box of ROI? In this case use Mexico
N = 33.
S = 14.
E = -118.
W = -86.

tile_list=np.genfromtxt(datadir+'tile_list.txt',dtype='string')
n_tiles = len(tile_list)
for i in range(0,n_tiles):
    if io.is_GeoTIFF_within_bbox(datadir+tile_list[i],N,S,W,E):
        print tile_list[i]

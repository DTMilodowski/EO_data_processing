# This function plots up a series of maps illustrating  forest loss through time,
# as modelled based on global forest watch annual forest loss distributed temporally
# based on the "seasonality" of FORMA 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np

import sys
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)

# Some basic info for maps
N = 33.
S = 14.
E = -86.
W = -118.

dY = 0.125
dX = 0.125
lat = np.arange(S,N,dY)#+dY/2. # shifting to cell centre
lon = np.arange(W,E,dX)#+dX/2. # shifting to cell centre

lon_grid,lat_grid = np.meshgrid(lon,lat)

# Load in GFW
GFW = np.load('GFW_monthly.npz')['arr_0']
N_months=GFW.shape[0]
cum_GFW = np.cumsum(GFW,axis=0)

ForestLoss=np.ma.masked_where(GFW<=0,GFW)

# Now make the plots
for i in range(0,N_months):
    fig = plt.figure(1, facecolor='White',figsize=[5,8])
    ax1a= plt.subplot2grid((2,1),(0,0)) 
    ax1a.set_title('Monthly forest loss')
    m1a = Basemap(projection='aea', lat_0=(N+S)/2., lon_0=(E+W)/2., llcrnrlat=S, urcrnrlat=N,llcrnrlon=W, urcrnrlon=E, resolution='i')
    m1a.ax = ax1a
    x,y = m1a(lon_grid,lat_grid)
    im1 = m1a.pcolormesh(x,y,GFW[i,:,:],vmin=0.0,vmax=0.004, rasterized=True, edgecolor='0.6', linewidth=0)
    cbar = m1a.colorbar(im1)
    cbar.solids.set_edgecolor("face")
    cbar.set_ticks([0,0.002,0.004])

    m1a.drawcountries(color='0.6',linewidth=1)
    m1a.drawcoastlines(color='0.5',linewidth=1)


    ax1b= plt.subplot2grid((2,1),(1,0)) 
    ax1b.set_title('Cumulative forest loss')
    m1b = Basemap(projection='aea', lat_0=(N+S)/2., lon_0=(E+W)/2., llcrnrlat=S, urcrnrlat=N,llcrnrlon=W, urcrnrlon=E, resolution='i')
    m1b.ax = ax1b
    x,y = m1b(lon_grid,lat_grid)
    im2=m1b.pcolormesh(x,y,cum_GFW[i,:,:],vmin=0.0,vmax=1, rasterized=True, edgecolor='0.6', linewidth=0)
    m1b.drawcountries(color='0.6',linewidth=1)
    m1b.drawcoastlines(color='0.5',linewidth=1)
    cbar = m1b.colorbar(im2)
    cbar.solids.set_edgecolor("face")
    cbar.set_ticks([0,0.5,1])

    plt.savefig('plot_EO_data/ForestLossMonthly_tstep'+str(i).zfill(3)+'.png')

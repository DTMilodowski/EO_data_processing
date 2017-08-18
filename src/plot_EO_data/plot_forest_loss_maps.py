# This function plots up a series of maps illustrating  forest loss through time,
# as modelled based on global forest watch annual forest loss distributed temporally
# based on the "seasonality" of FORMA 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np

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
regrid = np.load('GFW_annual.npz')['arr_0']
N_years=regrid.shape[0]
cum_regrid = np.cumsum(regrid,axis=0)

ForestLoss=np.ma.masked_where(regrid<=0,regrid)

# Now make the plots
for i in range(0,N_years):
    fig = plt.figure(1, facecolor='White',figsize=[5,8])
    ax1a= plt.subplot2grid((2,1),(0,0)) 
    ax1a.set_title('Monthly forest loss')
    m1a = Basemap(projection='aea', lat_0=(N+S)/2., lon_0=(E+W)/2., llcrnrlat=S, urcrnrlat=N,llcrnrlon=W, urcrnrlon=E, resolution='i')
    m1a.ax = ax1a
    x,y = m1a(lon_grid,lat_grid)
    m1a.pcolormesh(x,y,ForestLoss[i,:,:],vmin=0.0,vmax=0.2, cmap='jet', rasterized=True, edgecolor='0.6', linewidth=0)
    #cbar = m1a.colorbar()
    #cbar.solids.set_edgecolor("face")
    #cbar.set_ticks([0,0.1,0.2])

    m1a.drawcountries(color='0.6',linewidth=1)
    m1a.drawcoastlines(color='0.5',linewidth=1)


    ax1b= plt.subplot2grid((2,1),(1,0)) 
    ax1b.set_title('Cumulative forest loss')
    m1b = Basemap(projection='aea', lat_0=(N+S)/2., lon_0=(E+W)/2., llcrnrlat=S, urcrnrlat=N,llcrnrlon=W, urcrnrlon=E, resolution='i')
    m1b.ax = ax1b
    x,y = m1b(lon_grid,lat_grid)
    m1b.pcolormesh(x,y,cum_regrid[i,:,:],vmin=0.0,vmax=0.5, cmap='jet', rasterized=True, edgecolor='0.6', linewidth=0)
    m1b.drawcountries(color='0.6',linewidth=1)
    m1b.drawcoastlines(color='0.5',linewidth=1)

    plt.show()

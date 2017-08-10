# This function plots up a series of maps illustrating  forest loss through time,
# as modelled based on global forest watch annual forest loss distributed temporally
# based on the "seasonality" of FORMA 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import numpy as np

# Some basic info for maps
N = 33.
S = 14.
E = -86.
W = -118.


# Now make the plots
fig = plt.figure(1, facecolor='White',figsize=[5,8])
ax1a= plt.subplot2grid((2,1),(0,0)) 
ax1a.set_title('Monthly forest loss')
m1a = Basemap(projection='aea', lat_0=(N+S)/2., lon_0=(E+W)/2., llcrnrlat=S, urcrnrlat=N,llcrnrlon=W, urcrnrlon=E, resolution='i')
m1a.ax = ax1a
m1a.drawcountries(color='0.6',linewidth=1)
m1a.drawcoastlines(color='0.5',linewidth=1)


ax1b= plt.subplot2grid((2,1),(1,0)) 
ax1b.set_title('Cumulative forest loss')
m1b = Basemap(projection='aea', lat_0=(N+S)/2., lon_0=(E+W)/2., llcrnrlat=S, urcrnrlat=N,llcrnrlon=W, urcrnrlon=E, resolution='i')
m1b.ax = ax1b
m1b.drawcountries(color='0.6',linewidth=1)
m1b.drawcoastlines(color='0.5',linewidth=1)

plt.show()

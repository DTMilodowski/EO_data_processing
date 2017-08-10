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
fig = plt.figure(1, facecolor='White',figsize=[8,5])
ax1a= plt.subplot2grid((1,2),(0,0)) 

m1a = Basemap(projection='aea', lat_0=(N+S)/2., lon_0=(E+W)/2., llcrnrlat=S, urcrnrlat=N,llcrnrlon=W, urcrnrlon=E, resolution='i')
m1a.ax = ax1a


m1a.drawcountries(color='0.5',linewidth=1)
m1a.drawcoastlines(color='0.2',linewidth=1)
#m1a.fillcontinents()
plt.show()

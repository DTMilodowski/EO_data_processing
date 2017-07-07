# A bunch of functions that are useful for processing EO data
import numpy as np
#------------------------------------------------------------------------------
# Is point inside polygon? Ray-trace algorithm which utilises the numpy arrays
# so that a list of x and y coordinates can be processed in one call and only
# points inside polygon are returned alongside the indices in case required for
# future referencing.
def points_in_poly(x,y,poly):
    n = len(poly)
    inside=np.zeros(x.size,dtype=bool)
    xints=np.zeros(x.size)

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y=poly[i % n]
        if p1y!=p2y:
            xints[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)] = (y[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
        if p1x==p2x:
            inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)] = np.invert(inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x)],axis=0)])
        else:
            inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x),x<=xints],axis=0)] = np.invert(inside[np.all([y>min(p1y,p2y), y<=max(p1y,p2y), x<=max(p1x,p2x),x<=xints],axis=0)])
        p1x,p1y = p2x,p2y

    return x[inside],y[inside], inside



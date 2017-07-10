# import required python libraries 
import numpy as np

# a set of functions to calculate the area of a WGS84 pixel analytically
# for WGS84 a = 6378137 metres; b = 6356752.3142 metres
def calculate_area_of_ellipsoidal_slice(lat,a=6378137., b=6356752.3142):

    #convert lat to radians
    f = lat/360.*2*np.pi

    # calculating surface area from equator to latitude
    e = np.sqrt(1 - (b/a)**2)
    zm = 1 - e*np.sin(f)
    zp = 1 + e*np.sin(f)
    area = np.pi * b**2 * (np.log(zp/zm) / (2*e) + np.sin(f) / (zp*zm))
    return area

# returns area in square metres
def calculate_WGS84_pixel_area(lat1,lat2,long1,long2):

    # distance between east and west boundaries = fraction of whole circle
    q = np.abs(long1-long2)/360.

    # get surface areas from equator to each latitude
    a1 = calculate_area_of_ellipsoidal_slice(lat1)
    a2 = calculate_area_of_ellipsoidal_slice(lat2)

    # difference of these areas multiplied by fraction of 360 degrees taken up by
    # segment bounded by longitude values gives area
    area = q * (np.max([a1,a2]) - np.min([a1,a2]))

    return area

# this function produces an array of cell areas according to a specified range of latitudes and longitudes
# default units are sq. metres.  To change, specify area_scalar; e.g. for hactares -> area_scalar = 1./10.**4
# Also have the option to distunguish between providing lat and long for cell midpoint (default), or lower
# left corner
def calculate_cell_area(lat,long,area_scalar=1.,cell_centred=True):

    dx = long[1]-long[0]
    dy = lat[1]-lat[0]

    # shift lat and long so that they refer to cell boundaries if necessary
    if cell_centred == True:
        long=long-dx/2.
        lat = lat-dy/2.

    longs,lats = np.meshgrid(long,lat)
    rows,cols = longs.shape


    cell_area = np.zeros((rows,cols))
    for rr in range (0,rows):
        for cc in range(0,cols):
            lat1 = lats[rr,cc]
            lat2 = lat1+dy
            long1 = longs[rr,cc]
            long2 = long1+dx
            cell_area[rr,cc] = calculate_WGS84_pixel_area(lat1,lat2,long1,long2)

    # convert cell_area from sq. metres to desired units using scalar 
    cell_area*=area_scalar
    return cell_area


def calculate_cell_area_array(lat,long,area_scalar = 1.,cell_centred=True):
    dx = long[1]-long[0]
    dy = lat[1]-lat[0]

    # shift lat and long so that they refer to cell boundaries if necessary
    if cell_centred == True:
        long=long-dx/2.
        lat = lat-dy/2.

    longs,lats = np.meshgrid(long,lat)
    rows,cols = longs.shape

    q = np.abs(dx)/360.
    a1 =  calculate_area_of_ellipsoidal_slice(lats)
    a2 =  calculate_area_of_ellipsoidal_slice(lats+dy)

    cell_area = q * (np.max([a1,a2],axis = 0) - np.min([a1,a2],axis = 0))
    cell_area*=area_scalar
    return cell_area

# new function since cell area not a function of longitude, so only need a column
def calculate_cell_area_column(lat,dx,area_scalar = 1.,cell_centred=True):

    dy = lat[1]-lat[0]

    # shift lat and long so that they refer to cell boundaries if necessary
    if cell_centred == True:
        lat = lat-dy/2.

    q = np.abs(dx)/360.
    a1 =  calculate_area_of_ellipsoidal_slice(lat)
    a2 =  calculate_area_of_ellipsoidal_slice(lat+dy)

    cell_area = q * (np.max([a1,a2],axis = 0) - np.min([a1,a2],axis = 0))
    cell_area*=area_scalar
    return cell_area

#------------------------------------------------------------------------------
# Clip array to given bbox extent
# - arguments:
#   - input array
#   - input geotransformation
#   - bbox N limit
#   - bbox S limit
#   - bbox W limit
#   - bbox E limit
# - returns:
#   - clipped array
#   - geotransform for new clipped array
def clip_array_to_bbox(array,geoTrans,N,S,W,E):

    rows,cols = array.shape
    latitude = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3],geoTrans[5])
    longitude =  np.arange(geoTrans[0],cols*geoTrans[1]+geoTrans[0],geoTrans[1])
    d_lat = np.abs(geoTrans[5])
    d_long = np.abs(geoTrans[1])

    lat_keep = np.empty((latitude.size,1))
    long_keep = np.empty((1,longitude.size))

    # pad by half resolution for cases where pixels don't align perfectly with boundaries.
    lat_keep[:,0] = np.all((latitude <= N + d_lat/2., latitude >= S - d_lat/2.),axis=0)
    long_keep[0,:] = np.all((longitude <= E + d_long/2., longitude >= W - d_long/2.),axis=0)

    test = np.dot(lat_keep,long_keep).astype(bool)

    clip_array = array[test].reshape(lat_keep.sum(),long_keep.sum())
    
    geoTrans_u = []
    if geoTrans[5] > 0: # +ve y resolution indicates that raster origin specified as lower left corner
        geoTrans_u = [np.min(longitude[long_keep[0,:].astype(bool)]), geoTrans[1], geoTrans[2], np.min(latitude[lat_keep[:,0].astype(bool)]), geoTrans[4], geoTrans[5]]
    else:              # -ve y resolution indicates that raster origin specified as upper left corner
        geoTrans_u = [np.min(longitude[long_keep[0,:].astype(bool)]),  geoTrans[1], geoTrans[2], np.max(latitude[lat_keep[:,0].astype(bool)]), geoTrans[4], geoTrans[5]]
    
    return clip_array, geoTrans_u


#------------------------------------------------------------------------------
# Regrid array using nearest neighbour
# - arguments:
#   - target latitude (1D array)
#   - target longitude (1D array)
#   - initial geotransformation info
#   - initial array
#   - optional host array, for exxample if you are analysing tiled data sequentially
#   - resampling method; only option currently is 'sum' (default)
# - returns:
#   - regridded array
def regrid_array_nearest(target_lat, target_long, geoTrans, array, regrid = np.array([]), method = 'sum'):

    target_rows = target_lat.size
    target_cols = target_lon.size
    
    if regrid.size == 0:
        regrid = np.zeros((target_rows,target_cols))*np.nan
    
    init_rows,init_cols = array.shape
    closest_lat=np.zeros(init_rows).astype("int")
    closest_long=np.zeros(init_cols).astype("int")
    init_lat = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3]-0.000001*geoTrans[5],geoTrans[5])+geoTrans[5]/2. # shifting to cell centre
    init_long = np.arange(geoTrans[0],cols*geoTrans[1]+geoTrans[0]-0.000001*geoTrans[1],geoTrans[1])+geoTrans[1]/2. # shifting to cell centre

    for ii,val in enumerate(init_lat):
        closest_lat[ii]=np.argsort(np.abs(val-target_lat))[0]
    for jj,val in enumerate(init_long):
        closest_long[jj]=np.argsort(np.abs(val-target_long))[0]
    
    for ii in range(0,target_rows):
        lat_mask = closest_lat==ii
        if lat_mask.sum() > 0:
            for jj in range(0,target_cols):
                long_mask = closest_long==jj
                if long_mask.sum() > 0:
                    # avoid overwriting existing data from other tiles
                    if np.isnan(regrid[ii,jj]):
                        regrid[ii,jj] = np.sum(array[np.ix_(lat_mask,long_mask)])
                    else:
                        regrid[ii,jj] += np.sum(array[np.ix_(lat_mask,long_mask)])

    return regrid

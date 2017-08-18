# Functions to read and write EO data layers

import numpy as np
from matplotlib import pyplot as plt
from osgeo import gdal
import os
import osr
import sys
import datetime as dt
from netCDF4 import Dataset

# import my own libraries
import auxilliary_functions as aux

###############################################################################
# GeoTIFF functions
#------------------------------------------------------------------------------
# Function to read geoTIFF header and find whether it contains data within
# region of interest
def is_GeoTIFF_within_bbox(File,T,B,L,R):

    target_bbox = np.asarray([[L,T],[R,T],[R,B],[L,B]])

    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    try:
        ds = gdal.Open(File)
    except RuntimeError, e:
        print 'unable to open ' + File
        print e
        sys.exit(1)
    geoTrans = ds.GetGeoTransform()

    Xdim = ds.RasterXSize + 1 # since we need to traverse final pixel
    Ydim = ds.RasterYSize + 1 # since we need to traverse final pixel

    test  = False
    X0 = geoTrans[0]
    dX = geoTrans[1]
    Y0 = geoTrans[3]
    dY = geoTrans[5]

    BB11 = [X0 + Xdim*dX, Y0 + Ydim*dY]
    BB10 = [X0 + Xdim*dX, Y0]
    BB00 = [X0, Y0]
    BB01 = [X0 , Y0 + Ydim*dY]

    GeoTIFF_bbox = np.asarray([BB11,BB10,BB00,BB01])

    x,y,inside = aux.points_in_poly(GeoTIFF_bbox[:,0],GeoTIFF_bbox[:,1],target_bbox)
    if inside.sum()>0:
        test = True
    return test

#------------------------------------------------------------------------------
# Function to load a GeoTIFF band plus georeferencing information.  
# Only loads one band, which is band 1 by default
def load_GeoTIFF_band_and_georeferencing(File,band_number=1):
    
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    try:
        ds = gdal.Open(File)
    except RuntimeError, e:
        print 'unable to open ' + File
        print e
        sys.exit(1)
        
    source_band = ds.GetRasterBand(band_number)
    if source_band is None:
        print "BAND MISSING"
        sys.exit(1)  

    array = np.array(ds.GetRasterBand(band_number).ReadAsArray(),dtype=np.float64)
    geoTrans = ds.GetGeoTransform()
    coord_sys = ds.GetProjectionRef()

    return array, geoTrans, coord_sys

#------------------------------------------------------------------------------
# Function to write an array to a geoTIFF
def write_array_to_GeoTiff(array,geoTrans,OUTFILE_prefix,EPSG_CODE='4326',north_up=True):
    NBands = 1
    NRows = 0
    NCols = 0

    if north_up:
        # for north_up array, need the n-s resolution (element 5) to be negative
        if geoTrans[5]>0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2: 
            print 'array has less than two dimensions! Unable to write to raster'
            sys.exit(1)  
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print 'array has too many dimensions! Unable to write to raster'
            sys.exit(1)  

    else:
        # for north_up array, need the n-s resolution (element 5) to be positive
        if geoTrans[5]<0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
        # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
        if len(array.shape) < 2: 
            print 'array has less than two dimensions! Unable to write to raster'
            sys.exit(1)  
        elif len(array.shape) == 2:
            (NRows,NCols) = array.shape
            array = np.flipud(array)
        elif len(array.shape) == 3:
            (NRows,NCols,NBands) = array.shape
            for i in range(0,NBands):
                array[:,:,i] = np.flipud(array[:,:,i])
        else:
            print 'array has too many dimensions! Unable to write to raster'
            sys.exit(1)  
    
    # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
    if len(array.shape) < 2: 
        print 'array has less than two dimensions! Unable to write to raster'
        sys.exit(1)  
    elif len(array.shape) == 2:
        (NRows,NCols) = array.shape
        array = np.flipud(array)
    elif len(array.shape) == 3:
        (NRows,NCols,NBands) = array.shape
        for i in range(0,NBands):
            array[:,:,i] = np.flipud(array[:,:,i])
    else:
        print 'array has too many dimensions! Unable to write to raster'
        sys.exit(1)  
    
    # Write GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( OUTFILE_prefix+'.tif', NCols, NRows, NBands, gdal.GDT_Float32 )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    dataset.GetRasterBand(1).SetNoDataValue( -9999 )
    dataset.GetRasterBand(1).WriteArray( array )
    dataset = None
    return 0


###############################################################################
# netCDF functions
#------------------------------------------------------------------------------
def load_TRMM_NetCDF(NetCDF_file):

    ds = Dataset(NetCDF_file)
    pptn = np.transpose(ds.variables['precipitation'])
    lat =  np.transpose(ds.variables['lat'])
    lon =  np.transpose(ds.variables['lon'])

    return pptn, lat, lon


###############################################################################
# FORMA functions
# FORMA is rather annoyingly in a ascii file, which is a pain, but we'll have
# to deal with this for now
# This section is modified from some original code written by J-F Exbrayat
#------------------------------------------------------------------------------
# This function loads FORMA and resamples to a target resolution (the original
# resolution is 0.00426666667 deg
# It returns a multi-timestep grid with the fraction of the grid cell disturbed
# in the previous 16 days
def grid_FORMA(FORMAfile, target_resolution):

    original_resolution=0.00426666667

    lat=np.arange(90-target_resolution/2.,-90,-target_resolution)
    lon=np.arange(-180+target_resolution/2.,180,target_resolution)

    #open the ASCII dataset
    forma_data=file(FORMAfile,'r')
    #read the header line and the first data line
    line=forma_data.readline()
    line=forma_data.readline()

    #create the temporary lists
    ascii_dates=[]
    #iterate until the end of file is reached
    while line != '':
        date=line.split(',')[-1].strip()
        if date not in ascii_dates:
            ascii_dates.append(date)
        line=forma_data.readline()

    #sort dates
    ascii_dates.sort()

    #create the array to store data
    degradation=np.zeros([len(ascii_dates),len(lat),len(lon)])

    #go back to first line of file
    forma_data.seek(0)
    #read the header line and the first data line
    line=forma_data.readline()
    line=forma_data.readline()
    #iterate until the end of file is reached
    counter = 0
    while line != '':

        #find out the dates
        line=line.split(',')
        date=line[-1].strip()

        #get time step id
        idstep=ascii_dates.index(date)

        ptlat=float(line[0]);ptlon=float(line[1])
        idlat=np.argsort(np.abs(ptlat-lat))[0]
        idlon=np.argsort(np.abs(ptlon-lon))[0]
        
        #weight as function of cos for latitude
        weight=np.cos(np.radians(ptlat))/np.cos(np.radians(lat[idlat]))      
        degradation[idstep,idlat,idlon]=degradation[idstep,idlat,idlon]+(original_resolution/target_resolution)*weight*(original_resolution/target_resolution)
        
        line=forma_data.readline();counter+=1
        if counter % 100000 == 0:
            print counter

    forma_data.close()

    #create the time dimension
    timesteps=np.empty(len(ascii_dates),'i')

    for dd, date in enumerate(ascii_dates):

        day=np.array(date.split('-'),dtype='i')
        day=dt.datetime(day[0],day[1],day[2])
 
        if dd==0:
            refday=day        

        timesteps[dd]=(day-refday).days

    return timesteps, lat, lon, degradation


def grid_FORMA_monthly(FORMAfile,target_resolution,N,S,E,W,start_date = '2006-01-01', end_date='2015-01-01'):
    
    timesteps, lat, lon, degrad_i = grid_FORMA(FORMAfile,target_resolution)

    lat_mask = np.all((lat<N,lat>=S),axis=0)
    lon_mask = np.all((lon<E,lon>=W),axis=0)
    
    FORMAstart = np.datetime64('2005-12-19')
    
    # create array of dates for original FORMA dataset
    dates = FORMAstart+timesteps.astype('timedelta64[D]')
    
    # Clip the area to ROI
    degrad = np.zeros((dates.size,lon_mask.sum(),lat_mask.sum()))
    for dd in range(0,dates.size):
        degrad_iter = degrad_i[dd,:,:]
        degrad[dd,:,:] = degrad_iter[np.ix_(lat_mask,lon_mask)] # need to loop through for this as it isn't working at the moment
    degrad[degrad<0.] = 0.
    degrad_i=None
    degrad_iter=None


    months = np.unique(dates.astype('datetime64[M]')[1:])# want to skip the first  of the layers as this is funky
    N_months = months.size 
    monthly_degrad=np.zeros([N_months,degrad.shape[1],degrad.shape[2]])

    #start from degradation detected on or after 17/1/2006
    refmonth=dates[1].astype('datetime64[M]')-dates[1].astype('datetime64[Y]')
    monthid=0
    for dd,cday in enumerate(dates[1:]):
        print cday, dd, dates.size, monthid, N_months
        #get month of current time step
        currentmonth=dates[dd+1].astype('datetime64[M]')-dates[dd+1].astype('datetime64[Y]')
        # same month, just add data        
        if currentmonth == refmonth:
            print "same month", cday
            monthly_degrad[monthid] = monthly_degrad[monthid]+degrad[dd+1]
        #different month, distribute equally
        else:
            #find first day of month
            begin_month = cday.astype('datetime64[M]').astype('datetime64[D]')
            ndays = (cday-dates[dd+1-1]).astype('int') # number of days 
            nday_in_mth = (cday-begin_month).astype('int') #cday.day-1 #number of days in same month as date of record, this way we assume degradation did not occur on day of record
            nday_prev_mth = ndays-nday_in_mth 
            print "different month", cday, nday_in_mth, nday_prev_mth, ndays

            #assign portion of degradation to previous month
            monthly_degrad[monthid] = monthly_degrad[monthid]+degrad[dd+1]*nday_prev_mth/float(ndays)
            monthid = monthid +1
            #assign portion of degradation to currentmonth
            monthly_degrad[monthid] = monthly_degrad[monthid]+degrad[dd+1]*nday_in_mth/float(ndays)
            refmonth = currentmonth
            
    mask = np.all((months>=np.datetime64(start_date),months<np.datetime64(end_date)),axis=0)
    months = months[mask]
    monthly_degrad = monthly_degrad[mask,:,:]
    return months, lat_mask, lon_mask, monthly_degrad

    

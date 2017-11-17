#===============================================================================
# load_GFED.py
#-------------------------------------------------------------------------------
# @author D. T. Milodowski, November 2017
# This is a set of functions to load in GFED4 burned area data into numpy arrays 
# These arrays have three dimensions: lat, long, time, with a monthly timestep
#
# Data references:
# - Giglio, L., J.T. Randerson, and G.R. van der Werf, (2013), J. Geophys. Res. 
#   Biogeosci., 118, 317328, doi:10.1002/jgrg.20042.
# - Randerson, J.T., Y. Chen, G.R. van derWerf, B.M. Rogers, and D.C. Morton
#   (2012), J. Geophys. Res., 117, G04012, doi:10.1029/2012JG002128.
# - van der Werf, G.R., Randerson, J.T., Giglio, L., et al., (2017), Earth Syst.. 
#   Sci. Data, 9, 697-720, https://doi.org/10.5194/essd-9-697-2017.
#===============================================================================
# import standard libraries
import numpy as np

# import netcdf libraries
from netCDF4 import Dataset

# Function to load GFED4 monthly data. There are two potential variables here:
# - burned area expressed as fraction of a pixel: BurnedFraction
# - burned area expressed as the area in square metres: BurnedArea
# Files collate annual data, with timestep indicated by the day of the year.
# Need to specify time period of interest and the N,S,E,W extents for the area
# of interest.
def load_GFED4_monthly(path2files,variable,start_month,start_year,end_month,end_year,N,S,E,W):

    # first of all obtain the start and end date of the time series
    start_date = np.datetime64('%04i-%02i' % (start_year,start_month))
    end_date = (np.datetime64('%04i-%02i' % (end_year,end_month))+np.timedelta64(1,'M'))
    dates = np.arange(start_date,end_date)
    n_dates = dates.size

    # Load one tile to dimensions of clipped array, and the array mask
    NetCDF_file = '%s/GFED4_%04i.nc' % (path2files,start_year)
    ds = Dataset(NetCDF_file)
    lat = np.asarray(ds.variables['latitude'])
    lon = np.asarray(ds.variables['longitude'])
    
    lat_mask = np.all((lat<=N,lat>=S),axis=0)
    lon_mask = np.all((lon<=E,lon>=W),axis=0)
    n_lat = lat_mask.sum()
    n_lon = lon_mask.sum()
    ds=None
    
    # Loop through the netcdf files, retrieving the data and putting into time series.
    year = np.arange(start_year,end_year + 1)
    month = np.arange(12)+1
    n_years = year.size
    i_mm = 0
    GFED_sample = np.zeros((n_dates,n_lat,n_lon))
    for yy in range(0,n_years):
        
        # get the number of months needed for this year
        n_months = 12
        start_mm = 1
        end_mm = 12
        if yy == 0:
            n_months = 12 - start_month + 1
            start_mm = start_month
        if year[yy] == end_year:
            n_months = n_months - (12-end_month)
            end_mm = end_month
        
        NetCDF_file = '%s/GFED4_%04i.nc' % (path2files,year[yy])
        print NetCDF_file

        ds = Dataset(NetCDF_file)
        # get area of interest
        month_mask = np.all((month>=start_mm,month<=end_mm),axis=0)
        print month_mask.shape, lat_mask.shape, lon_mask.shape
        array_mask = np.ix_(month_mask,lat_mask,lon_mask)

        GFED_sample[i_mm:i_mm+n_months] = np.asarray(ds.variables[variable])[array_mask]
        i_mm += n_months

    return dates, GFED_sample

#%% [markdown]
# These functions are designed to download ERA5 data per month. In other words, the function will split the required duration into 1-month intervals. The idea behind this module is to download files in two steps:
# 
# 1. Request the data to the server, which will queue all requests. All requests are stored in `out_path/ERA5props.csv`, but no download is made yet. In that way, a single python session doesn't need to be permanently running!
# 2. Manually check what requests have been completed and download those that have. Again, this updates `ERA5props.csv` so only new files are downloaded
# 
# ### References:
# 
# - https://github.com/ecmwf/cdsapi/issues/2
# - https://github.com/ecmwf/cdsapi/blob/master/examples/example-era5-update.py

# %%
import cdsapi
import calendar
from pathlib import Path
import pandas as pd
import webbrowser
import numpy as np
import xarray as xr
from windrose import WindroseAxes
import matplotlib.cm as cm

__all__ = [
    'plevelList', 'dayList', 'hourList', 'set_area', 'make_date_vec', 'submitERA', 
    'checkERA', 'queue', 'loadNc', 'ncTodf', 'getPoint', 'vec2wind', 'getLevel', 'findNearest'
]

plevelList = ['1', '2', '3','5', '7', '10','20', '30', '50','70', '100', '125','150', '175', '200','225', '250', '300','350', '400', '450','500', '550', '600','650', '700', '750','775', '800', '825','850', '875', '900','925', '950', '975','1000']
dayList = ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30','31']
hourList = ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']

# %%

def set_area(lat, lon, nx=1, ny=1, res=0.25, single_point=False):
    """
    Create a bounding box for a grid centered on the closest ERA5 grid point to (lat, lon)
    and expands it by (nx, ny) and resolution res (e.g., nx=1 and ny=1 will return a 3x3 grid).

    Parameters:
        lat (float): Latitude of the target point.
        lon (float): Longitude of the target point.
        nx (int): Number of grid points by which to expand the center point in the longitudinal direction (default 1).
        ny (int): Number of grid points by which to expand the center point in the latitudinal direction (default 1).
        res (float): Grid resolution in degrees (default 0.25 for ERA5).
        single_point (bool): If True, return the bounding box as a single point (default False).

    Returns:
        list: Bounding box coordinates [North, West, South, East].
    """
    
    if single_point:
        return [lat, lon, lat, lon]
    else:
        # Snap to the closest grid point
        center_lat = np.round(lat / res) * res
        center_lon = np.round(lon / res) * res

        # Calculate half-extent based on grid size
        lat_extent = (ny * res) / 2
        lon_extent = (nx * res) / 2

        # Define the bounding box
        north = center_lat + lat_extent
        south = center_lat - lat_extent
        west = center_lon - lon_extent
        east = center_lon + lon_extent

        return [north, west, south, east]

def make_date_vec(year_start, year_end, month_start, month_end, day_start=None, day_end=None):
    """
    Generate a list of date vectors based on the specified start and end years, months, and days.

    Args:
        year_start (int): The starting year.
        year_end (int): The ending year.
        month_start (int): The starting month.
        month_end (int): The ending month.
        day_start (int, optional): The starting day. Defaults to None.
        day_end (int, optional): The ending day. Defaults to None.

    Returns:
        list: A list of date vectors, where each vector contains [day, month, year] if day_start and day_end are specified,
              or [month, year] if day_start and day_end are None.
    """
    
    stor = []
    # Loop through years
    for iY in range(year_start, year_end+1):
        # If first year, start with `month_start`
        if (iY==year_start) & (month_start>1):
            iMs=month_start
        # For following years, start in January
        else:
            iMs=1
            
        # Similarly, take `month_end` for the last year and December otherwise
        if (iY==year_end) & (month_end<13):
            iMe=month_end
        else:
            iMe=12
        
        # Loop through months 
        for iM in range(iMs, iMe+1):
            
            # If days are not specified, the vector will contain [month, year]
            if day_start == None:
                stor.append([iM, iY])
            
            # Otherwise, loop through the days
            else:
                
                # If first month start with `day_start`, otherwise take 1
                if iM == month_start:
                    iDs = day_start
                else:
                    iDs = 1
                
                # Same
                if iM == month_end:
                    iDe = day_end 
                else:
                    lastday = calendar.monthrange(year_start,month_start)
                    iDe = lastday[1]
                
                # Loop through days
                for iD in range(iDs, iDe+1):
                    stor.append([iD, iM, iY])
    return stor

def submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, day_start=1, day_end=31, format='netcdf', dt='month'):
    """
    Submits a CDS API request for each month or day of the requested time interval.

    Parameters:
    - out_path (str): The path to the output folder where the downloaded data will be saved.
    - year_start (int): The starting year of the time interval.
    - year_end (int): The ending year of the time interval.
    - month_start (int): The starting month of the time interval.
    - month_end (int): The ending month of the time interval.
    - dataset (str): The name of the dataset to retrieve from the CDS API.
    - rDict (dict): A dictionary containing the parameters for the CDS API request.
    - day_start (int, optional): The starting day of the time interval. Defaults to 1.
    - day_end (int, optional): The ending day of the time interval. Defaults to 31.
    - format (str, optional): The format of the downloaded data. Defaults to 'netcdf'.
    - dt (str, optional): The time interval of the request. Can be 'month' or 'day'. Defaults to 'month'.

    Returns:
    None
    """
    
    # Make sure the folder exists
    Path(out_path).mkdir(parents=True, exist_ok=True)

    # Create the date vector   
    if dt == 'month':
        date_vec = make_date_vec(year_start, year_end, month_start, month_end)
    else:
        date_vec = make_date_vec(year_start, year_end, month_start, month_end, day_start=day_start, day_end=day_end)

    # Check that last day actually matches what is specified in year/month 2022-03-18
    day_end = day_end if day_end==calendar.monthrange(year_end,month_end)[1] else calendar.monthrange(year_end,month_end)[1]
    
    # Setup storage dataframe
    df = pd.DataFrame()
     
    # Loop through combinations of months/years and submit queries
    count  = 1
    for iDate in date_vec:
        
        if dt == 'month':
            month, year = iDate[0], iDate[1]
            print(f'Accessing {dataset} data for {month}/{year}')
            
            # Computes days 2022-03-18
            dS = day_start if month==month_start else 1
            dE = day_end if month==month_end else calendar.monthrange(year,month)[1]
            
            # 2022-03-18 using /to/ was returning an error when not using the `-complete` dataset required for WRF. This seems to fix it
            if dataset == 'reanalysis-era5-complete':
                rDict['date'] = f'{year}{month:02d}{dS:02d}/to/{year}{month:02d}{dE:02d}'
            else:
                rDict['date'] = f'{year}{month:02d}{dS:02d}/{year}{month:02d}{dE:02d}'
        else:
            day, month, year = iDate[0], iDate[1], iDate[2]
            print(f'Accessing {dataset} data for {day}/{month}/{year}')
            
            # 2022-11-02 Trying Alex's approach
            rDict['date'] = f'{year}{month:02d}{day:02d}'
            
        if format == 'netcdf':
            rDict['format'] = 'netcdf'
        
        print(rDict)
        # Start request
        c = cdsapi.Client(wait_until_complete=False, delete=False)
        r = c.retrieve(dataset, rDict)

        # Update the storage
        if dt == 'month':
            df = pd.concat([df, pd.DataFrame({'r': r.reply['request_id'], 'month': month, 'year': year, 'format': format, 'completed': False}, index=[count])])
        else:
            df = pd.concat([pd.DataFrame({'r': r.reply['request_id'], 'day': day, 'month': month, 'year': year, 'format': format, 'completed': False}, index=[count])])
            
        count += 1
    
    # Save storage to file
    df.to_csv(f'{out_path}/ERA5props.csv')
    
def checkERA(out_path, prepend=False):
    """
    Checks the status of ERA5 data files and downloads any pending files.

    Args:
        out_path (str): The path to the output directory.
        prepend (bool): Whether to prepend the file name with file number.

    Returns:
        None
    """
    df = pd.read_csv(f'{out_path}/ERA5props.csv', index_col=0)
    
    # Check if all files have already been downloaded
    if df.completed.sum() == df.shape[0]:
        print('All files have already been downloaded and I refuse to work overtime...')
    else:
        cnt = 0
        for i in range(0, df.shape[0]):
            # Extract information from the dataframe
            count = df.iloc[i].name
            year = df.iloc[i].year
            month = df.iloc[i].month
            format = df.iloc[i]['format']
            
            if 'day' in df.columns:
                day = df.iloc[i].day
       
            if format == 'netcdf':
                fileExt = 'nc'
            else: 
                fileExt = 'grib'
                
            # Generate the output file path
            if prepend:
                if 'day' in df.columns:
                    outfile = "%s/%05d_%04d_%02d_%02d.%s"%(out_path, count, year, month, day, fileExt)
                else:
                    outfile = "%s/%05d_%04d_%02d.%s"%(out_path, count, year, month, fileExt)
            else:
                if 'day' in df.columns:
                    outfile = "%s/%04d_%02d_%02d.%s"%(out_path, year ,month,day, fileExt)
                else:
                    outfile = "%s/%04d_%02d.%s"%(out_path, year ,month, fileExt)
            
            # Check if the file is already completed
            if not df.iloc[i].completed:
                reply = dict(request_id=df.iloc[i]['r']) # request_id from above
                new_client = cdsapi.Client(wait_until_complete=False, delete=False)
                # result = cdsapi.api.Result(new_client, reply)
                result = new_client.client.get_remote(df.iloc[i]['r']) 
                result.update()
                reply = result.reply
                print(result)
                if reply['state'] == 'completed':
                    cnt += 1
                    df.loc[df.iloc[[i]].index, 'completed'] = True
                    df.to_csv(f'{out_path}/ERA5props.csv') 
                    try:
                        result.download(outfile)
                    except:
                        print(f'Error downloading {outfile}')
        
        # Check if all files have been completed
        if df.completed.sum() == df.shape[0]:
            print('Woohooo - all done!')
        else:
            print(f'{cnt} new files completed, total = {df.completed.sum() }/{df.shape[0]}. For more info, call `queue`')
        
def queue():
    webbrowser.open('https://cds.climate.copernicus.eu/requests?tab=all', new=0, autoraise=True)

def loadNc(out_path, save_path=None, time_dim='valid_time', process=True, drop_uv=True, dtype=np.int16):
    """
    Load NetCDF files from the specified `out_path` and optionally save the loaded data to a new NetCDF file at `save_path`.

    Parameters:
    - out_path (str): The path to the NetCDF files to be loaded.
    - save_path (str, optional): The path to save the loaded data as a new NetCDF file. If not provided, the data will not be saved.
    - time_dim (str, optional): The name of the time dimension in the NetCDF files. Defaults to 'valid_time'.
    - process (bool, optional): Whether to process the loaded data (e.g., compute wind speed and direction). Defaults to True.
    - drop_uv (bool, optional): Whether to drop the original u and v wind components after computing wind speed and direction. Defaults to True.
    
    Returns:
    - xarray.Dataset: The loaded wind data.

    """
    def calculate_scale_offset(data, dtype=np.int16):
        """Calculate optimal scale_factor and add_offset for a given data array."""
        data_min = np.nanmin(data)
        data_max = np.nanmax(data)
        
        # Handle cases with no range
        if data_max == data_min:
            return 1.0, 0.0
        
        # Set add_offset as the midpoint
        add_offset = (data_max + data_min) / 2
        
        # Calculate scale_factor to fit data in the integer range
        int_info = np.iinfo(dtype)
        scale_factor = (data_max - data_min) / (int_info.max - int_info.min)
        
        return scale_factor, add_offset
    
    # Load data
    ds = xr.open_mfdataset(f'{out_path}*.nc',decode_cf=False)
    # Sort by time
    ds = ds.sortby(time_dim)
    # Make sure time is named 'time'
    if time_dim != 'time':
        ds = ds.rename({time_dim: 'time'})
    
    # Compute altitude, speed and direction
    if process:
        print('- Processing data')
        # Compute wind speed and direction
        ds = vec2wind(ds, drop_uv=drop_uv)
    
    # Clean some variables
    if 'number' in ds.data_vars:
            ds = ds.drop('number')
    if 'expver' in ds.data_vars:
            ds = ds.drop('expver')
            
    # Set compression options
    encoding = {}
    for var in ds.data_vars:
        data = ds[var].values
        
        # Compute scale and offset
        scale_factor, add_offset = calculate_scale_offset(data, dtype)
        
        # Store the encoding options
        encoding[var] = {
            "dtype": dtype.__name__,
            "scale_factor": scale_factor,
            "add_offset": add_offset,
            "zlib": True,
            "complevel": 4,
        }

    # Save if requested
    if save_path is not None:
        print(f'- Saving data to {save_path}')
        ds.load().to_netcdf(save_path, encoding=encoding)
    
    return ds
    
    
def ncTodf(ds):
    """
    Converts an xarray dataset to a pandas dataframe and reorders the indices.

    Parameters:
    ds (xarray.Dataset): The xarray dataset to be converted.

    Returns:
    pandas.DataFrame: The converted pandas dataframe with reordered indices.
    """
    
    # Convert xarray to pandas and reorder index
    df = ds.to_dataframe()

    # Reorder the indices to query with loc. Need to get rid of `level` in case of ERA5 Land
    if 'level' in df.index.names:
        levels_var = ['latitude', 'longitude', 'level', 'time']
    else:
        levels_var = ['latitude', 'longitude', 'time']
    
    # In case ERA5 Land, rename `u10` to `u`
    if 'u10' in df.columns:
        df = df.rename({'u10':'u','v10':'v'},axis=1)
    df = df.round(2)
        
    return df.reorder_levels(levels_var)
    
def getPoint(df, lat, lon):

    # Make sure lon is between 0 and 360
    if lon<0: lon = 360 + lon
    
    # Find closest point
    lonVal = findNearest(df.index.get_level_values('longitude').values, lon)
    latVal = findNearest(df.index.get_level_values('latitude').values, lat)

    # Extract dataframe for closest point
    if 'level' in df.index.names:
        # print(df.head())
        return df.loc[latVal, lonVal, :, :]
    else:
        # print(df.head())
        return df.loc[latVal, lonVal, :]
    
def vec2wind(ds, drop_uv=True):
    """
    Converts vector components (u and v) to wind speed and direction.

    Parameters:
    - ds (pandas.DataFrame or xarray.Dataset): The input data containing vector components.
    - drop_uv (bool, optional): Whether to drop the u and v columns after conversion. Default is True.

    Returns:
    - ds (pandas.DataFrame or xarray.Dataset): The input data with added columns for wind speed and direction.
    """
    
    # Get elevation in m asl
    ds['altitude'] = ds['z']/9.80665 
    # Compute wind direction
    ds['wind_direction'] = (np.arctan2(ds['v'], ds['u']) * 180 / np.pi) % 360
    # Compute wind speed
    ds['wind_speed'] = np.sqrt(ds['u']**2 + ds['v']**2)

    # Drop u and v if requested
    if drop_uv:
        if isinstance(ds, xr.core.dataset.Dataset):
            ds = ds.drop(['u', 'v'])
        else:
            ds = ds.drop(['u', 'v'], axis=1)
        
    return ds

def getLevel(df, alt):

    # Compute mean elevation and join that back to df
    meanZ = df['z'].groupby('level').mean().rename('z_mean')
    df = df.join(meanZ)

    # Find the closest elevation
    altVal = findNearest(df.z_mean.values, alt)
    dfAlt = df[df.z_mean == altVal]
    
    return dfAlt


def findNearest(array, num):
    # nearest_idx = np.where(abs(array-num)==abs(array-num).min())[0] # If you want the index of the element of array (array) nearest to the the given number (num)
    # nearest_val = array[abs(array-num)==abs(array-num).min()] # If you directly want the element of array (array) nearest to the given number (num)
    val = np.unique(array[abs(array-num)==abs(array-num).min()]) # If you directly want the element of array (array) nearest to the given number (num)
    return val[0]

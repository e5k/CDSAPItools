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
import numpy as np
import xarray
from windrose import WindroseAxes
import matplotlib.cm as cm


plevelList = ['1', '2', '3','5', '7', '10','20', '30', '50','70', '100', '125','150', '175', '200','225', '250', '300','350', '400', '450','500', '550', '600','650', '700', '750','775', '800', '825','850', '875', '900','925', '950', '975','1000'],
dayList = ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30','31']
hourList = ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']

# %%
def make_date_vec(year_start, year_end, month_start, month_end, day_start=None, day_end=None):
    """ Create a list of list containing all permutations of [month, year]"""
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
    """ Submits a CDS API request for each month of the requested time interval. 
    
     Args:
        out_path (str): Path where files will be downloaded
        year_start (int): Start year
        year_end (int): End year
        month_start (int): Start month
        month_end (int): End month
        dataset (str): A valid CDS dataset, e.g., `reanalysis-era5-pressure-levels` or `reanalysis-era5-land` (see https://confluence.ecmwf.int/display/CKB/The+family+of+ERA5+datasets)
        rDict (dict): Dictionary to be sent to the API
        day_start (int): First day to retrieve (default: 1)
        day_end (int): Last day to retrieve (default: 31)
        format (str): Accepts 'grib' and 'netcdf'
        dt (str): Defines the temporal granularity, accepts 'month' and 'day'
        
    Returns:
        pd.DataFrame: A file used to later check the status and download requested files is saved under `out_path/ERA5props.csv`

    """
    
    # Make sure the folder exists
    Path(out_path).mkdir(parents=True, exist_ok=True)

    # Create the date vector   
    if dt == 'month':
        date_vec = make_date_vec(year_start, year_end, month_start, month_end)
    else:
        date_vec = make_date_vec(year_start, year_end, month_start, month_end, day_start=day_start, day_end=day_end)
        
    # Setup area
    if rDict['area'][1]<0: rDict['area'][1]=360+rDict['area'][1]
    if rDict['area'][3]<0: rDict['area'][3]=360+rDict['area'][3]

    # Setup storage dataframe
    df = pd.DataFrame()
     
    # Loop through combinations of months/years and submit queries
    count  = 1
    for iDate in date_vec:
        
        if dt == 'month':
            month, year = iDate[0], iDate[1]
            print(f'Accessing ERA5 data for {month}/{year}')
            rDict['date'] = f'{year}{month:02d}01/to/{year}{month:02d}31'
        else:
            day, month, year = iDate[0], iDate[1], iDate[2]
            print(f'Accessing ERA5 data for {day}/{month}/{year}')
            rDict['date'] = f'{year}{month:02d}{day:02d}/to/{year}{month:02d}{day:02d}'
        
        if format == 'netcdf':
            rDict['format'] = 'netcdf'
        
        # Start request
        c = cdsapi.Client(wait_until_complete=False, delete=False)
        r = c.retrieve(dataset, rDict)

        # Update the storage
        if dt == 'month':
            df = df.append(pd.DataFrame({'r': r.reply['request_id'], 'month': month, 'year': year, 'format': format, 'completed': False}, index=[count]))
        else:
            df = df.append(pd.DataFrame({'r': r.reply['request_id'], 'day': day, 'month': month, 'year': year, 'format': format, 'completed': False}, index=[count]))
            
        count += 1
    
    # Save storage to file
    df.to_csv(f'{out_path}/ERA5props.csv')
    
def checkERA(out_path, prepend):
    
    df = pd.read_csv(f'{out_path}/ERA5props.csv', index_col=0)
    
    if df.completed.sum() == df.shape[0]:
        print('All files have already been downloaded and I refuse to work overtime...')
    else:
    #     print(f'{df.completed.sum() }/{df.shape[0]} files downloaded...')
        
        cnt = 0
        for i in range(0, df.shape[0]):
            
            count = df.iloc[i].name
            year = df.iloc[i].year
            month = df.iloc[i].month
            format = df.iloc[i]['format']
            
            if format == 'netcdf':
                fileExt = '.nc'
            else: 
                fileExt = '.grb'
                
            if prepend:
                outfile = "%s/%05d_%04d_%02d.nc"%(out_path, count, year, month)
            else:
                outfile = "%s/%04d_%02d.nc"%(out_path, year ,month)
            
            if not df.iloc[i].completed:
                reply = dict(request_id=df.iloc[i]['r']) # request_id from above
                new_client = cdsapi.Client(wait_until_complete=False, delete=False)
                result = cdsapi.api.Result(new_client, reply)
                result.update()
                reply = result.reply
                print(result)
                if reply['state'] == 'completed':
                    cnt += 1
                    df.iloc[i]['completed'] = True
                    result.download(outfile)
        
        df.to_csv(f'{out_path}/ERA5props.csv')            
        
        if df.completed.sum() == df.shape[0]:
            print('Woohooo - all done!')
        else:
            print(f'{cnt} new files completed, total = {df.completed.sum() }/{df.shape[0]}. For more info, call `queue`')
        
def queue():
    webbrowser.open('https://cds.climate.copernicus.eu/cdsapp#!/yourrequests', new=0, autoraise=True)

def loadNc(path):
    
    # Load data
    ds = xarray.open_dataset(path)

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
        
    return df.reorder_levels(levels_var)
    
def getPoint(df, lat, lon):

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
    
def vec2wind(df):   
    # Compute velocity and direction from u and v
    df['speed'] = np.sqrt(np.power(df['u'],2) + np.power(df['v'],2))
    df['direction'] = np.degrees(np.arctan2(df['u'],df['v']))
    df.loc[df['direction']<0, 'direction'] = 360 + df.loc[df['direction']<0, 'direction']

    return df

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
# %%

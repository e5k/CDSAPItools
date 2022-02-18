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

# %%
def make_date_vec(year_start, year_end, month_start, month_end):
    """ Create a list of list containing all permutations of [month, year]"""
    stor = []
    for iY in range(year_start, year_end+1):
        if (iY==year_start) & (month_start>1):
            iMs=month_start
        else:
            iMs=1
        if (iY==year_end) & (month_end<13):
            if month_end == month_start:
                iMe = month_start+1
            else:
                iMe=month_end
        else:
            iMe=12
            
        for iM in range(iMs, iMe+1):
            stor.append([iM, iY])
            
    return stor
   
def submitIt(dataset, pressure_level, year, month, dayList, hourList, grid, variables, area):            
    
    lastday1=calendar.monthrange(year,month)
    lastday=lastday1[1]
    dayList = range(lastday+1)
    dayList = dayList[1:]
    dayList = [str(i) for i in dayList]

    bdate="%s%02d01"%(year,month)
    edate="%s%02d%s"%(year,month,lastday)

    print('Accessing ERA5 data from ', bdate,' to ',edate,' (YYYYMMDD)')


    c = cdsapi.Client(wait_until_complete=False, delete=False)
    r = c.retrieve(
        dataset, 
        {
            'variable'      : variables,
            'pressure_level': pressure_level,
            'product_type'  : 'reanalysis',
            'year'          : '%s'%(year),
            'month'         : '%s'%(month),
            'day'           : dayList,       
            'area'          : area, # North, West, South, East. Default: global
            'grid'          : grid, # Latitude/longitude grid: east-west (longitude) and north-south resolution (latitude). Default: 0.25 x 0.25
            'time'          : hourList,
            'format'        : 'netcdf' # Supported format: grib and netcdf. Default: grib
        })

    return r

def submitERA(out_path, year_start, year_end, month_start, month_end, north, south, west, east, dataset, prepend, pressure_level, dayList, hourList, grid, variables):
    
    # Make sure the folder exists
    Path(out_path).mkdir(parents=True, exist_ok=True)

    # Create the date vector   
    date_vec = make_date_vec(year_start, year_end, month_start, month_end)
        
    # Setup area
    if west<0: west=360+west
    if east<0: east=360+east
    area = [north, west, south, east]

    df = pd.DataFrame()
     
    # Download
    count  = 1
    for iDate in date_vec:
        month, year = iDate[0], iDate[1]

        

        r = submitIt(dataset, pressure_level, year, month, dayList, hourList, grid, variables, area)
        
        df = df.append(pd.DataFrame({'r': r.reply['request_id'], 'month': month, 'year': year, 'completed': False}, index=[count]))
        
        count += 1
        
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
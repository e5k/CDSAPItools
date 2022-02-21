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

plevelList = ['1', '2', '3','5', '7', '10','20', '30', '50','70', '100', '125','150', '175', '200','225', '250', '300','350', '400', '450','500', '550', '600','650', '700', '750','775', '800', '825','850', '875', '900','925', '950', '975','1000'],
dayList = ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30','31']
hourList = ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']

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

def submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict):
    """ Submits a CDS API request for each month of the requested time interval. 
    
     Args:
        out_path (str): Path where files will be downloaded
        year_start (int): Start year
        year_end (int): End year
        month_start (int): Start month
        month_end (int): End month
        dataset (str): A valid CDS dataset, e.g., `reanalysis-era5-pressure-levels` or `reanalysis-era5-land` (see https://confluence.ecmwf.int/display/CKB/The+family+of+ERA5+datasets)
        rDict (dict): Dictionary to be sent to the API

    Returns:
        pd.DataFrame: A file used to later check the status and download requested files is saved under `out_path/ERA5props.csv`

    """
    
    # Make sure the folder exists
    Path(out_path).mkdir(parents=True, exist_ok=True)

    # Create the date vector   
    date_vec = make_date_vec(year_start, year_end, month_start, month_end)
        
    # Setup area
    if rDict['area'][1]<0: rDict['area'][1]=360+rDict['area'][1]
    if rDict['area'][3]<0: rDict['area'][3]=360+rDict['area'][3]
    # area = [north, west, south, east]

    # Setup storage dataframe
    df = pd.DataFrame()
     
    # Loop through combinations of months/years and submit queries
    count  = 1
    for iDate in date_vec:
        month, year = iDate[0], iDate[1]

        # Setup day list
        # lastday1=calendar.monthrange(year,month)
        # lastday=lastday1[1]
        # dayList = range(lastday+1)
        # dayList = dayList[1:]
        # dayList = [str(i) for i in dayList]

        # bdate="%s%02d01"%(year,month)
        # edate="%s%02d%s"%(year,month,lastday)

        print(f'Accessing ERA5 data from {month}/{year}')

        # Update dictionary
        rDict['year'] = f'{year}'
        rDict['month'] = f'{month}'
        rDict['format'] = 'netcdf'
        
        # Start request
        c = cdsapi.Client(wait_until_complete=False, delete=False)
        r = c.retrieve(dataset, rDict)

        # Submit job
        # r = submitIt(dataset, pressure_level, year, month, dayList, hourList, grid, variables, area)

        # Update the storage
        df = df.append(pd.DataFrame({'r': r.reply['request_id'], 'month': month, 'year': year, 'completed': False}, index=[count]))
        
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
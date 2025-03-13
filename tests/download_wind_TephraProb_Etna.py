# %%
import cdsapitools.main as cds
import schedule

#%% Etna
# Set spatio-temporal properties of dataset:
# Time of dataset
year_start  = 2014
year_end    = 2023
month_start = 1
month_end   = 12

# Define area [north, west, south, east]
area = cds.set_area(37.750320, 14.993779, single_point=True)    

# Set some global parameters
out_path = '/Users/seb/Documents/WORK/Students/UCLouvain/Master/Zoe Saintrain/wind-etna-2014-2023'
hourList = [x+':00' for x in ['00', '06', '12', '18']]

# Choose the dataset
dataset = 'reanalysis-era5-pressure-levels'

# Choose the variables as a function of the selected dataset. 
variableList = ['geopotential', 'u_component_of_wind', 'v_component_of_wind']
grid = [0.25, 0.25]
dataset = 'reanalysis-era5-pressure-levels'

# Adapt the request dictionary to reflect the chosen dataset. For instance, ERA5 land doesn't need the `pressure_level` field. If you are not sure of what this dictionary is supposed to look like, visit the reference for any dataset provided above, click on the `Download data` tab, setup a mock download and on the `Show API request` green button at the end of the page.
# Do not add the `year` or `month` field as it will be dynamically added later
# Note: default values are available as `plevelList`, `dayList`, `hourList`
rDict = {
    "product_type": "reanalysis",
    "format": "netcdf",
    "variable": cds.variableList,
    "pressure_level": cds.plevelList,
    "time": hourList,
    "area": area,
    "grid": grid,
}

#%% Submit the request
cds.submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict)

#%% Check the status
cds.checkERA(out_path)

#%% Open the website to see the queue
cds.queue()

#%% If you want to use the scheduler to check updates daily
schedule.every().day.at("14:30").do(cds.checkERA(out_path))
while True:
    schedule.run_pending()
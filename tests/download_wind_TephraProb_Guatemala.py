# %%
import cdsapitools.main as cds
import schedule

#%% Guatemala
# Set spatio-temporal properties of dataset:
# Time of dataset
year_start  = 1990
year_end    = 1990
month_start = 1
month_end   = 2

# Define area [north, west, south, east]
area = [15.5, -92.5, 13.5, -89.5]

# Set some global parameters
out_path = '/var/services/homes/seb/ERA/Guatemala1900-2000'
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
    "variable": variableList,
    "pressure_level": cds.plevelList,
    "time": hourList,
    "area": area,
    "grid": grid,
}

#%% Submit the request
cds.submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict)

#%% Check the status
# cds.checkERA(out_path)

#%% Open the website to see the queue
# cds.queue()

#%% If you want to use the scheduler to check updates daily
schedule.every().day.at("14:30").do(cds.checkERA(out_path))
while True:
    schedule.run_pending()
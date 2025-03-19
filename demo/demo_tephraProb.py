#%% [markdown]
# # TephraProb demo
# This notebook illustrates how to use `cdsapitools` to download ERA5 data and use it to run TephraProb simulations. It uses a case study at Etna volcano and retrieves 4-daily ERA5 data for the period 2014-2023. As the cds api limits the number of points per request, `cdsapitools` splits the request per months.

# Make sure you have setup the CDS API as described [here](https://cds.climate.copernicus.eu/how-to-api)
#%%
import cdsapitools.main as cds
import schedule

#%% [markdown]
# ## Setup the configuration
# We need to setup a configuration dictionary `rDict` to be passed to the `cds.submitERA` function. This dictionary is used to define the request to the CDS API. If you are not sure of what this dictionary is supposed to look like, visit the reference for any ERA5, click on the `Download data` tab, setup a mock download and on the `Show API request` green button at the end of the page.

#%%
# Set the output path
out_path = '/Users/wind-etna-2014-2023'

# Set the temporal extent of the data to be retrieved. `cdsapitools` retrieves one file per year/month, which will be concatenated later.
year_start  = 2014
year_end    = 2023
month_start = 1
month_end   = 12

# Set the hours to retrieve. You can use `cds.hourList` to retrieve all hours or set it manually.
hourList = [x+':00' for x in ['00', '06', '12', '18']]

# Define area in the shape [north, west, south, east]. The `set_area` function is a wrapper to estimate areas or points. See the documentation for more details. In this case, we retrieve only one point, which means that wind interpolation should be turned off in TephraProb.
area = cds.set_area(37.750320, 14.993779, single_point=True)

# Set the vertical grid. You can use `cds.plevelList` to retrieve all levels or set it manually.
pLevels = cds.plevelList

# Set the grid size
grid = [0.25, 0.25]

# Choose the dataset
dataset = 'reanalysis-era5-pressure-levels'

# Choose the variables as a function of the selected dataset. Adapt the requested variables to reflect the chosen dataset. For instance, ERA5 land doesn't need the `pressure_level` field. 
variableList = ['geopotential', 'u_component_of_wind', 'v_component_of_wind']

# Setup the main config dictionary. Do not add the `year` or `month` field as it will be dynamically added later
# Note: default values are available as `plevelList`, `dayList`, `hourList`
rDict = {
    "product_type": "reanalysis",
    "format": "netcdf",
    "variable": variableList,
    "pressure_level": pLevels,
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
# schedule.every().day.at("14:30").do(cds.checkERA(out_path))
# while True:
#     schedule.run_pending()

#%% If some jobs have failed, you can resubmit them
cds.submitFailed(out_path, dataset, rDict)

#%% Concatenate the result files and save them to a new file. For use in TephraProb, use process=False and drop_uv=False
ds = cds.loadNc(out_path+'/', '/Users/seb/Documents/WORK/Students/UCLouvain/Master/Zoe Saintrain/ERA5_Etna_2014-2023.nc', process=False, drop_uv=False)


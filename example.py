# %% [markdown]
# # CDSAPItools
# 
# These functions are designed to download ERA5 data per month. In other words, the function will split the required duration into 1-month intervals. The idea behind this module is to prevent having to keep a Python instance running during the entire duration of the queue process. This is implemented in two steps:
# 
# 1. Request the data to the server, which will queue all requests. All requests are stored in `out_path/ERA5props.csv`, but no download is made yet. 
# 2. Manually check what requests have been completed and download those that have. Again, this updates `ERA5props.csv` so only new files are downloaded
# 
# 
# ## Requirements
# 
# `CDSAPItools` require the `pandas`, `cdsapi` and `windrose` modules:
# 
# ```bash
# pip install cdsapi pandas windrose
# ````
#  
# ## References:
# 
# ### CDS background
# 
# - Homepage: https://cds.climate.copernicus.eu/#!/home
# - Datasets: https://cds.climate.copernicus.eu/cdsapp#!/search?type=dataset
# 
# 
# ### Method implementation:
# 
# - https://github.com/ecmwf/cdsapi/issues/2
# - https://github.com/ecmwf/cdsapi/blob/master/examples/example-era5-update.py
# 
# 
# ## To do:
# 
# - [x] 2022-02-17: Still need to adapt the queries in `submitIt` to match different datasets - for now only `reanalysis-era5-pressure-levels`
#   - 2022-02-21: Partially fixed, now accepts `reanalysis-era5-pressure-levels` and `reanalysis-era5-land`
# - [ ] 2022-02-18: Find why counter is not updating

# %%
# If not running from the same folder, uncomment these two lines and adapt the path
import os   
os.chdir('/Users/seb/Documents/Codes/CDSAPItools')
from CDSAPItools import *
import pandas as pd


#%% [markdown]
# ## Set global parameters

# %%
# Set spatio-temporal properties of dataset:
# Time of dataset
year_start  = 2021
year_end    = 2021
month_start = 9
month_end   = 12

# Define area [north, west, south, east]. W and E are in degrees E, so:
#   if west<0: west=360+west
#   if east<0: east=360+east
area = [28.7, 342.05, 28.5, 342.25]

# Set some global parameters
# Output folder, i.e. replace by your project name
out_path    = '/Users/seb/Documents/WORK/Projects/La Palma/Wind/'

# Prepend file number to output file name (useful for TephraProb)
prepend = False

# Select dataset
# Define the dataset: https://confluence.ecmwf.int/display/CKB/The+family+of+ERA5+datasets
# e.g.:
# - 'reanalysis-era5-pressure-levels': https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview
# - 'reanalysis-era5-land': https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

#%% [markdown] 
# ## Example 1: Pressure levels

# %%
dataset = 'reanalysis-era5-pressure-levels'
out_path = '/Users/seb/Documents/WORK/Projects/La Palma/Wind/PL'

# Choose the variables as a function of the selected dataset. 
# Note: default values are available as `plevelList`, `dayList`, `hourList`
plevelList = ['200','225', '250', '300','350', '400', '450','500', '550', '600','650', '700', '750','775', '800', '825','850', '875', '900','925', '950', '975','1000']
dayList = ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30','31']
hourList = ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']
variableList = ['geopotential', 'u_component_of_wind', 'v_component_of_wind']
grid = [0.25, 0.25]

# Adapt the request dictionary to reflect the chosen dataset. For instance, ERA5 land doesn't need the `pressure_level` field. If you are not sure of what this dictionary is supposed to look like, visit the reference for any dataset provided above, click on the `Download data` tab, setup a mock download and on the `Show API request` green button at the end of the page.
# Do not add the `year` or `month` field as it will be dynamically added later
rDict = {
    "product_type": "reanalysis",
    "format": "netcdf",
    "variable": variableList,
    "pressure_level": plevelList,
    "day": dayList,
    "time": hourList,
    "area": area,
    "grid": grid,
}

#%% [markdown] 
# ## Example 2: Surface

# %%
dataset = 'reanalysis-era5-land'
out_path = '/Users/seb/Documents/WORK/Projects/La Palma/Wind/SF'

# Choose the variables as a function of the selected dataset. 
# Note: default values are available as `plevelList`, `dayList`, `hourList`
dayList = ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30','31']
hourList = ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']
variableList = ['10m_u_component_of_wind', '10m_v_component_of_wind']
grid = [0.1, 0.1]

# Adapt the request dictionary to reflect the chosen dataset. For instance, ERA5 land doesn't need the `pressure_level` field. If you are not sure of what this dictionary is supposed to look like, visit the reference for any dataset provided above, click on the `Download data` tab, setup a mock download and on the `Show API request` green button at the end of the page.
# Do not add the `year` or `month` field as it will be dynamically added later

rDict = {
    "format": "netcdf",
    "variable": variableList,
    "day": dayList,
    "time": hourList,
    "area": area,
    "grid": grid,
}


# %%
# submitERA(out_path, year_start, year_end, month_start, month_end, north, south, west, east, dataset, prepend, plevelList, dayList, hourList, grid, variableList)
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict)

# %%
checkERA(out_path, prepend)

# %% Concatenate all monthly data into a time-continuous dataset
# https://neetinayak.medium.com/combine-many-netcdf-files-into-a-single-file-with-python-469ba476fc14

import netCDF4
import numpy as np
import xarray

ds = xarray.open_mfdataset(f'{out_path}/*.nc',combine = 'by_coords', concat_dim="time")

# Save as a new netCDF file
ds.to_netcdf(f'{out_path}/{year_start}{month_start:02d}-{year_end}{month_end:02d}.nc')

#%% Extract and plot the wind profile for a given point
# 

def findNearest(array, num):
    # nearest_idx = np.where(abs(array-num)==abs(array-num).min())[0] # If you want the index of the element of array (array) nearest to the the given number (num)
    # nearest_val = array[abs(array-num)==abs(array-num).min()] # If you directly want the element of array (array) nearest to the given number (num)
    return np.unique(array[abs(array-num)==abs(array-num).min()]) # If you directly want the element of array (array) nearest to the given number (num)

lat = 28.4
lon = -18



# 
if lon<0: lon=360+lon

# Convert xarray to pandas and reorder index
df = ds.to_dataframe()
df = df.reorder_levels(['latitude', 'longitude', 'level', 'time'])

# Compute velocity and direction from u and v
df['velocity'] = np.sqrt(np.power(df['u'],2) + np.power(df['v'],2))
df['direction'] = np.degrees(np.arctan2(df['u'],df['v']))
df.loc[df['direction']<0, 'direction'] = 360 + df.loc[df['direction']<0, 'direction']

# Find closest point
df_lon = df.index.get_level_values('longitude').values
df_lat = df.index.get_level_values('latitude').values
lonVal = findNearest(df_lon, lon)
latVal = findNearest(df_lat, lat)

# Extract dataframe for closest point
dfLoc = df.loc[latVal, lonVal, :, :]

# Compute mean elevation and join that back to df
meanZ = dfLoc['z'].groupby('level').mean().rename('z_mean')
dfLoc = dfLoc.join(meanZ)
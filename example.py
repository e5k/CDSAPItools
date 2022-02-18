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
# `CDSAPItools` require the `pandas` and `cdsapi` modules:
# 
# ```bash
# pip install cdsapi pandas
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
# - [ ] 2022-02-17: Still need to adapt the queries in `submitIt` to match different datasets - for now only `reanalysis-era5-pressure-levels`
# - [ ] 2022-02-18: Find why counter is not updating

# %%
# If not running from the same folder, uncomment these two lines and adapt the path
# import os   
# os.chdir('/Users/seb/Documents/Codes/CDSAPItools')
from CDSAPItools import *
import pandas as pd

# %%
# Time of dataset
year_start  = 2021
year_end    = 2021
month_start = 9
month_end   = 12

# Area
north       = 28.6+.1
south       = 28.6-.1
west        = -17.85-.1
east        = -17.85+.1

#%% Set parameters
pressure_level = ['200','225', '250', '300','350', '400', '450','500', '550', '600','650', '700', '750','775', '800', '825','850', '875', '900','925', '950', '975','1000']
# pressure_level = ['1', '2', '3','5', '7', '10','20', '30', '50','70', '100', '125','150', '175', '200','225', '250', '300','350', '400', '450','500', '550', '600','650', '700', '750','775', '800', '825','850', '875', '900','925', '950', '975','1000'],
dayList = ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30','31']
hourList = ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']
grid = [0.25, 0.25]
variables = ['geopotential', 'u_component_of_wind', 'v_component_of_wind']

#%% Choose dataset (default: 'reanalysis-era5-pressure-levels')
# 

# dataset = 'reanalysis-era5-land'
dataset = 'reanalysis-era5-pressure-levels'

# Output folder, i.e. replace by your project name
out_path    = '/Users/seb/Documents/WORK/Projects/La Palma/Wind/'

# Prepend file number to output file name 
prepend = False



# %%
submitERA(out_path, year_start, year_end, month_start, month_end, north, south, west, east, dataset, prepend, pressure_level, dayList, hourList, grid, variables)

# %%
checkERA(out_path, prepend)

# %% Concatenate all monthly data into a time-continuous dataset
# https://neetinayak.medium.com/combine-many-netcdf-files-into-a-single-file-with-python-469ba476fc14

import netCDF4
import numpy
import xarray

ds = xarray.open_mfdataset(f'{out_path}/*.nc',combine = 'by_coords', concat_dim="time")

# Save as a new netCDF file
ds.to_netcdf(f'{out_path}/{year_start}{month_start:02d}-{year_end}{month_end:02d}.nc')

# %%

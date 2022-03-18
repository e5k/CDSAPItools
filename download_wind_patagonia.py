# %% [markdown]


# %%
# If not running from the same folder, uncomment these two lines and adapt the path
import os   
os.chdir('/Users/seb/Documents/Codes/CDSAPItools')
from CDSAPItools import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray
from windrose import WindroseAxes
import matplotlib.cm as cm

# %% Concatenate all monthly data into a time-continuous dataset
# https://neetinayak.medium.com/combine-many-netcdf-files-into-a-single-file-with-python-469ba476fc14

# Read and concatenate all nc files
# ds = xarray.open_mfdataset(f'{out_path}/*.nc',combine = 'by_coords', concat_dim="time")

# # Save as a new netCDF file
# ds.to_netcdf(f'{out_path}/{year_start}{month_start:02d}-{year_end}{month_end:02d}.nc')

# %% Cordon Caulle / WRF Model level data

year_start  = 2006
year_end    = 2016
month_start = 5
month_end   = 7
area = [-39, -72.5, -43, -66.5] 
prepend = False
dataset = 'reanalysis-era5-complete'
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/CC2011/ML'

rDict = {
    'class':'ea',
    'area': area,
    'expver':'1',
    'levelist': '1/to/137',
    'levtype':'ml',
    'param': "129/130/131/132/133/152",
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid': "0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib')
checkERA(out_path, prepend)

# %% Cordon Caulle / WRF Surface data

paramList = "172/165/166/167/168/134/151/235/31/34/33/141/139/170/183/236/39/40/41/42"
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/CC2011/SFC'

rDict = {
    'class':'ea',
    'area':area,
    'expver':'1',
    'levtype':'sfc',
    'param': paramList, 
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid':"0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib')
checkERA(out_path, prepend)

# %% Chaiten/ WRF Model level data

year_start  = 2003
year_end    = 2013
month_start = 4
month_end   = 6
area = [-40.5, -73.0, -44.5, -67.5] 
prepend = False
dataset = 'reanalysis-era5-complete'
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/Chaiten2008/ML'

rDict = {
    'class':'ea',
    'area': area,
    'expver':'1',
    'levelist': '1/to/137',
    'levtype':'ml',
    'param': "129/130/131/132/133/152",
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid': "0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib')
checkERA(out_path, prepend)

# %% Chaiten/ WRF Surface data

paramList = "172/165/166/167/168/134/151/235/31/34/33/141/139/170/183/236/39/40/41/42"
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/Chaiten2008/SFC'

rDict = {
    'class':'ea',
    'area':area,
    'expver':'1',
    'levtype':'sfc',
    'param': paramList, 
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid':"0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib')
checkERA(out_path, prepend)

# %% Calbuco / WRF Model level data
   
year_start  = 2010
year_end    = 2020
month_start = 3
month_end   = 5
area = [-69, -41.5, -73, -38]
prepend = False
dataset = 'reanalysis-era5-complete'
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/Calbuco2015/ML'

rDict = {
    'class':'ea',
    'area': area,
    'expver':'1',
    'levelist': '1/to/137',
    'levtype':'ml',
    'param': "129/130/131/132/133/152",
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid': "0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib')
checkERA(out_path, prepend)

# %% Calbuco / WRF Surface data

paramList = "172/165/166/167/168/134/151/235/31/34/33/141/139/170/183/236/39/40/41/42"
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/Calbuco2015/SFC'

rDict = {
    'class':'ea',
    'area':area,
    'expver':'1',
    'levtype':'sfc',
    'param': paramList, 
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid':"0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib')
checkERA(out_path, prepend)



# %% Debug to add the daily granularity option

year_start  = 2013
year_end    = 2013
month_start = 6
month_end   = 6
day_start = 6
day_end = 7

area = [-39, -72.5, -43, -66.5] 
prepend = False
dataset = 'reanalysis-era5-complete'
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/CC2011/ML'

rDict = {
    'class':'ea',
    'area': area,
    'expver':'1',
    'levelist': '1/to/137',
    'levtype':'ml',
    'param': "129/130/131/132/133/152",
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid': "0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib',
          day_start=day_start, day_end=day_end, dt='day')
checkERA(out_path, prepend)

#%%

paramList = "172/165/166/167/168/134/151/235/31/34/33/141/139/170/183/236/39/40/41/42"
out_path = '/Users/seb/Documents/WORK/Projects/Vegetation/WRF/CC2011/SFC'

rDict = {
    'class':'ea',
    'area':area,
    'expver':'1',
    'levtype':'sfc',
    'param': paramList, 
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid':"0.25/0.25",
}
# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grib',
          day_start=day_start, day_end=day_end, dt='day')
checkERA(out_path, prepend)

# %%

# %%
# If not running from the same folder, uncomment these two lines and adapt the path
import os   
os.chdir('/Users/seb/Documents/Codes/CDSAPItools')
from CDSAPItools import *
import pandas as pd


# %% Tansen ###################################
# Set spatio-temporal properties of dataset:
# Time of dataset
year_start  = 2015
year_end    = 2016
month_start = 1
month_end   = 12

# Define area [north, west, south, east]
area = [-40, -74, -43, -71]

# Set some global parameters
# Output folder, i.e. replace by your project name
out_path    = '/Users/seb/Documents/WORK/Students/UNIGE Masters/Tansen Rahman/Wind1015-2016'

# Prepend file number to output file name (useful for TephraProb)
# prepend = True

# %% Guatemala ###################################
# Set spatio-temporal properties of dataset:
# Time of dataset
year_start  = 2010
year_end    = 2022
month_start = 1
month_end   = 12

# Define area [north, west, south, east]
area = [15.5, -92.5, 13.5, -89.5]

# Set some global parameters
# Output folder, i.e. replace by your project name
out_path    = '/Users/seb/Documents/CONTRACTS/IDMC/Datasets/Guatemala/ERA5'
hourList = [x+':00' for x in ['00', '06', '12', '18']]




# %% Etna ###################################
# Set spatio-temporal properties of dataset:
# Time of dataset
year_start  = 2014
year_end    = 2023
month_start = 1
month_end   = 12

# Define area [north, west, south, east]
area = set_area(37.750320, 14.993779, single_point=True)    

# Set some global parameters
out_path    = '/Users/seb/Documents/WORK/Students/UCLouvain/Master/Zoe Saintrain/wind-etna-2014-2023'
hourList = [x+':00' for x in ['00', '06', '12', '18']]

dataset = 'reanalysis-era5-pressure-levels'

# Choose the variables as a function of the selected dataset. 
# Note: default values are available as `plevelList`, `dayList`, `hourList`
variableList = ['geopotential', 'u_component_of_wind', 'v_component_of_wind']
grid = [0.25, 0.25]


# %%
dataset = 'reanalysis-era5-pressure-levels'

# Choose the variables as a function of the selected dataset. 
# Note: default values are available as `plevelList`, `dayList`, `hourList`
variableList = ['geopotential', 'u_component_of_wind', 'v_component_of_wind']
grid = [0.25, 0.25]


##% Geneva ###################################
# Set spatio-temporal properties of dataset:
# Time of dataset
year_start  = 2010
year_end    = 2017
month_start = 1
month_end   = 12

# Define area [north, west, south, east]
area = [46.25, 6, 46, 6.25]

# Set some global parameters
# Output folder, i.e. replace by your project name
out_path    = '/Users/seb/Documents/WORK/Students/UNIGE Masters/Tansen Rahman/Wind_GVA'
hourList = [x+':00' for x in ['00', '06', '12', '18']]

dataset = 'reanalysis-era5-pressure-levels'

# Choose the variables as a function of the selected dataset. 
# Note: default values are available as `plevelList`, `dayList`, `hourList`
variableList = ['geopotential', 'u_component_of_wind', 'v_component_of_wind', 'temperature', 'specific_humidity']


#%%
grid = [0.25, 0.25]

# Adapt the request dictionary to reflect the chosen dataset. For instance, ERA5 land doesn't need the `pressure_level` field. If you are not sure of what this dictionary is supposed to look like, visit the reference for any dataset provided above, click on the `Download data` tab, setup a mock download and on the `Show API request` green button at the end of the page.
# Do not add the `year` or `month` field as it will be dynamically added later
rDict = {
    "product_type": "reanalysis",
    "format": "netcdf",
    "variable": variableList,
    "pressure_level": plevelList,
    "time": hourList,
    "area": area,
    "grid": grid,
}


# %%
# submitERA(out_path, year_start, year_end, month_start, month_end, north, south, west, east, dataset, prepend, plevelList, dayList, hourList, grid, variableList)
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict)

# %%
checkERA(out_path)

#%%
# ds = loadNc(out_path+'/', '/Users/seb/Documents/WORK/Students/UNIGE Masters/Tansen Rahman/2015-2016.nc', process=True)
ds = loadNc(out_path+'/', '/Users/seb/Documents/CONTRACTS/IDMC/Datasets/Guatemala/ERA5_2010-2022.nc', process=True)

#%% Geneva
ds = loadNc(out_path+'/', '/Users/seb/Documents/WORK/Students/UNIGE Masters/Tansen Rahman/ERA5_Geneva_2010-2020.nc', process=True, drop_uv=False)

#%%

# # Correct directions
# if isinstance(df, xr.core.dataset.Dataset):
#     df = xr.where(df.direction<0, 360+df.direction, df.direction )
# else:
#     df.loc[df['direction']<0, 'direction'] = 360 + df.loc[df['direction']<0, 'direction']



# ds = xr.open_mfdataset(f'{out_path}/*.nc',combine = 'by_coords', concat_dim="time",decode_cf=False)
# ds = xr.open_mfdataset(f'{out_path}/*.nc',decode_cf=False)


# %% Concatenate all monthly data into a time-continuous dataset
# https://neetinayak.medium.com/combine-many-netcdf-files-into-a-single-file-with-python-469ba476fc14

import numpy as np
import xarray
from windrose import WindroseAxes
import matplotlib.cm as cm

# Read and concatenate all nc files
ds = xarray.open_mfdataset(f'{out_path}/*.nc',combine = 'by_coords', concat_dim="time")

# Save as a new netCDF file
ds.to_netcdf(f'{out_path}/{year_start}{month_start:02d}-{year_end}{month_end:02d}.nc')

#%% Extract and plot the wind profile for a given point
import matplotlib.pyplot as plt

lat = 28.4
lon = -17.7
alt = 1000

path = f'{out_path}/{year_start}{month_start:02d}-{year_end}{month_end:02d}.nc'

df = loadNc(path)
df = getPoint(df, lat, lon)
df = vec2wind(df)
df = getLevel(df, alt)

df1000 = getLevel(df, 1000).reset_index()
df2000 = getLevel(df, 3000).reset_index()

ht = pd.read_excel('/Users/seb/Documents/WORK/Projects/La Palma/Altura de la columna eruptiva.xlsx')

# %%
fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10,8))

ax1.plot(ht.time, ht.ht, 'k')

ax2.plot(df1000.time, df1000.speed)
ax2.plot(df2000.time, df2000.speed)

ax3.plot(df1000.time, df1000.direction)
ax3.plot(df2000.time, df2000.direction)


# %%
df1000 = df1000.reset_index().set_index('time')
df2000 = df2000.reset_index().set_index('time')
dfT1000 = df1000.resample('d').mean()
# dfT2000 = df2000.resample('d').mean()

# %%
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(20,15))
ax1.plot(ht.time, ht.ht, 'k')
ax2.quiver(dfT1000.reset_index().time, 1, dfT1000.u, dfT1000.v, width=.002)
ax3.plot(df1000.reset_index().time, df1000.speed)
ax4.plot(df1000.reset_index().time, df1000.direction)
# ax3.quiver(dfT2000.reset_index().time, 1, dfT2000.u, dfT2000.v, width=.002)

ax1.set_ylabel('Plume height (m)')
ax2.set_ylabel('Wind quiver')
ax3.set_ylabel('Wind speed (m/s)')
ax4.set_ylabel('Wind direction (deg)')




# %%
ds = xarray.open_dataset(f'{out_path}/{year_start}{month_start:02d}-{year_end}{month_end:02d}.nc')
if lon<0: lon=360+lon

# Convert xarray to pandas and reorder index
df = ds.to_dataframe()

# Reorder the indices to query with loc. Need to get rid of `level` in case of ERA5 Land
if 'level' in df.index.names:
    levels_var = ['latitude', 'longitude', 'level', 'time']
else:
    levels_var = ['latitude', 'longitude', 'time']
    
df = df.reorder_levels(levels_var)

# Compute velocity and direction from u and v
df['speed'] = np.sqrt(np.power(df['u'],2) + np.power(df['v'],2))
df['direction'] = np.degrees(np.arctan2(df['u'],df['v']))
df.loc[df['direction']<0, 'direction'] = 360 + df.loc[df['direction']<0, 'direction']

# Find closest point
lonVal = findNearest(df.index.get_level_values('longitude').values, lon)
latVal = findNearest(df.index.get_level_values('latitude').values, lat)

# Extract dataframe for closest point
dfPoint = df.loc[latVal, lonVal, :, :]
dfPoint = df.loc[(latVal, lonVal)]
dfPoint = df.loc[(latVal, lonVal),['u', 'v']]

# Compute mean elevation and join that back to df
meanZ = dfPoint['z'].groupby('level').mean().rename('z_mean')
dfPoint = dfPoint.join(meanZ)

# Find the closest elevation
altVal = findNearest(dfPoint.z_mean.values, alt)
dfAlt = dfPoint[dfPoint.z_mean == altVal]

# %%

# %%
ax = WindroseAxes.from_ax()
ax.bar(dfAlt.direction, dfAlt.speed, normed=True, opening=1, edgecolor='k', cmap=cm.viridis, bins=[2,5,10,15,20,30], blowto=True, linewidth=.25)
ax.set_legend()

# plot_windrose(dfAlt.reset_index().set_index(['time']), kind='bar', bins=np.arange(0.01,8,1), cmap=cm.viridis, lw=3, theta_labels=["E", "N-E", "N", "N-W", "W", "S-W", "S", "S-E"])
# %%









# %% Example for WRF Model level data
# https://dreambooker.site/2018/04/20/Initializing-the-WRF-model-with-ERA5/

year_start  = 2021
year_end    = 2021
month_start = 9
month_end   = 9

area = [28.7, 342.05, 28.5, 342.25]

# Prepend file number to output file name (useful for TephraProb)
prepend = False

dataset = 'reanalysis-era5-complete'
out_path = '/Users/seb/Documents/WORK/Projects/GEE/Wind/ML'

# Choose the variables as a function of the selected dataset. 
# Note: default values are available as `plevelList`, `dayList`, `hourList`
plevelList = "1/to/137"
# dayList = ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12','13', '14', '15','16', '17', '18','19', '20', '21','22', '23', '24','25', '26', '27','28', '29', '30','31']
# dayList = ['01', '02', '03']
hourList = ['00:00']
# hourList = ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']
paramList = "129"
# paramList = "129/130/131/132/133/152"
grid = "0.25/0.25"

rDict = {
    'class':'ea',
    'area': area,
    'expver':'1',
    'levelist': '1/to/137',
    'levtype':'ml',
    'param': paramList,
    'stream':'oper',
    'time': '00/to/23/by/1',
    'type':'an',
    'grid':grid,
}

# %% Example for WRF surface data
# https://dreambooker.site/2018/04/20/Initializing-the-WRF-model-with-ERA5/

submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, day_start=1, day_end=2, format='grib')
checkERA(out_path, prepend)
queue()

# %%
paramList = "172/165/166/167/168/134/151/235/31/34/33/141/139/170/183/236/39/40/41/42"
out_path = '/Users/seb/Documents/WORK/Projects/GEE/Wind/SFC'

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
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, day_start=1, day_end=2, format='grib')
checkERA(out_path, prepend)


# LEVELIST="1/to/137"
# PARAMSFC="172/165/166/167/168/134/151/235/31/34/33/141/139/170/183/236/39/40/41/42"
# PARAMML="129/130/131/132/133/152"


# CLASS="ea"
# DATASET="era5"
# EXPVER="1"
# GRID="0.25/0.25"
# STREAM="OPER"
# TIME="00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00"
# TYPE="an"


# mars <<EOFMARS
# RETRIEVE,
# CLASS = $CLASS,
# DATASET = $DATASET,
# TYPE = $TYPE,
# STREAM = $STREAM,
# EXPVER = $EXPVER,
# LEVTYPE = ml,
# LEVELIST = $LEVELIST,
# PARAM = $PARAMML,
# TIME = $TIME,
# AREA = $NORTH/$WEST/$SOUTH/$EAST,
# GRID = $GRID,
# DATE = $day,
# TARGET = "$target"
# %%

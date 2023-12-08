import netCDF4 as nc
import numpy as np

nc_data = nc.Dataset("data/S20RTS_dvs.nc", "r")
depths = nc_data.variables['depth'][:]
dvs = nc_data.variables['v'][:]
latitudes = nc_data.variables['latitude'][:]
longitudes = nc_data.variables['longitude'][:]

print(depths)

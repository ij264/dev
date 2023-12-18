import numpy as np
import pyshtools as pysh
import matplotlib.pyplot as plt
from tqdm import tqdm


# Read in data
data = np.loadtxt('res_depths.dat', delimiter=' ')

# Reshape data
# Column 0: longitude, Column 1: latitude, Column 2: residual depth, Column 3: error
data = data.reshape((int(data.shape[0]*4), int(data.shape[1]/4)))
# Sort data
data = data[np.argsort(data[:, 0])]


# gridl = pysh.SHGrid.from_array(data[:, 2].reshape((int(data.shape[0]/4), 4)))
longitudes = data[:, 0]
latitudes = data[:, 1]
sorted_latitudes = np.sort(latitudes)
res_depths = data[:, 2]
errors = data[:, 3]

griddh = np.zeros((len(latitudes), len(longitudes)))

for lon_index in range(len(longitudes)):
    # print(longitudes[lon_index], latitudes[lon_index], res_depths[lon_index])
    # find index of corresponding latitude
    lat_index = np.where(sorted_latitudes == latitudes[lon_index])
    griddh[lat_index, lon_index] = res_depths[lon_index]

clm = pysh.SHGrid.from_array(griddh).expand()

# Plot dynamic topography

longs, lats = np.meshgrid(longitudes, sorted_latitudes)

DT_grid = []
for lat in tqdm(lats):
    DT_grid.append(clm.expand(grid='DH', lat=lat, lon=longitudes))

DT_grid = np.array(DT_grid)

fig, ax = plt.subplots()
plt.imshow(DT_grid, cmap='RdBu_r', extent=[0, 360, -90, 90], origin='lower')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(extend='both', shrink=0.5)
plt.show()
# Calculate spherical harmonics
# clm = pysh.SHExpandDH(griddh, lmax_calc=20)
# print(data)
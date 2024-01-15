import numpy as np
import process
import pyshtools as pysh
import pygmt

# -----------------------------------PYGMT configuration-------------------------------
pygmt.config(COLOR_BACKGROUND="0.6375/0.6375/255", COLOR_FOREGROUND="255/0.6375/0.6375")
# -------------------------------------------------------------------------------------

# First we fit using degrees 1 to 40.
# We then want to plot degrees 1 to 3 of this, so we zero the other degrees.

l_max = 40
l_top = 40

clm_obs = pysh.SHCoeffs.from_file(fname='coefficients/degree_1_to_40.shcoef',
                                  format='shtools',
                                  csphase=-1,
                                  normalization='ortho',
                                  kind='real')

# Set the degrees greater than 3 to zero.
for l in range(l_top + 1, l_max + 1):
    for m in range(-l, l + 1):
        clm_obs.set_coeffs(values=0., ls=l, ms=m)

# Plot dynamic topography
latitudes = np.linspace(90, -90, 181)
longitudes = np.linspace(0, 360, 361)

clm_obs_arr = clm_obs.to_array()
clm_obs_grid = pysh.expand.MakeGrid2D(cilm=clm_obs_arr,
                                      interval=1,
                                      norm=4,
                                      csphase=-1)

clm_obs_grid = np.roll(clm_obs_grid, shift=180, axis=1)
np.savez_compressed(f'global_grid/degree_1_to_{l_top}_of_{l_max}', lats=latitudes, longs=longitudes, grid=clm_obs_grid)
xr_DT = process.np_to_xr(clm_obs_grid, lats=latitudes, lons=longitudes)

fig = pygmt.Figure()

print(pygmt.grdinfo(xr_DT))

with fig.subplot(
        nrows=1, ncols=1, figsize=("10c", "4c")):
    with fig.set_panel(panel=0):
        fig.grdimage(xr_DT,
                     cmap='coolwarm_DT.cpt',
                     shading=None,
                     projection='H12c',
                     frame=["a", f"t+tDynamic Topography"],
                     interpolation='b')
        fig.grdcontour(xr_DT,
                       interval=.500,
                       limit=[-2.000, 2.000],
                       frame=True,
                       projection='H12c')
        fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H12c',
                  region="d")

        fig.colorbar(cmap='coolwarm_DT.cpt', frame=["x+lElevation",
                                                 "y+lkm"])

print(f'Saving...')
fig.savefig(
    f'figures/maps/map_degree_1_to_{l_top}_of_{l_max}.png',
    dpi=300)

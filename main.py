import process
import numpy as np
import pyshtools as pysh
import matplotlib.pyplot as plt
import read_in
from tqdm import tqdm
import pygmt
import sys
import subprocess
import os

# Retrieve command-line arguments
if len(sys.argv) != 3:
    print("Usage: python main.py <viscosity> <kernel_type>")
    sys.exit(1)

valid_kernels = ['surftopo', 'geoid']

visc_src = sys.argv[1]
kernel_type = sys.argv[2]


assert isinstance(visc_src, str), "Kernel source must be a string."
assert kernel_type in valid_kernels, f"Kernel type must be one of {valid_kernels}."

# Your script logic using the arguments
print(f"Argument 1: {visc_src}")


if not os.path.exists(f'lyness/OUTPUT_{visc_src}/'):
# Replace 'your_command_here' with the command you want to run
    command = f'cd lyness && ./MAKE_KERNEL OUTPUT_{visc_src}/ VISC_INPUTS/const_visc.vis && cd ..'

    # Run the command in the shell
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(result.stdout)
        print("Kernels produced successfully")
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        print(e.stderr)

if not os.path.exists(f'figures/S20RTS/{visc_src}/{kernel_type}/'):
    os.makedirs(f'figures/S20RTS/{visc_src}/{kernel_type}/')

# Maximum spherical harmonic degree.
l_max = 20
# Read in data
depths, latitudes, longitudes, dvs = read_in.read_in_nc("data/S20RTS_dvs.nc")

# Note that longitude is from -180 to 180, so we need to shift it to 0 to 360.

# Convert to percentage.
dvs /= 100
# Convert to metres.
depths *= 1.e3

# Convert to shear wave velocity to density.
density_anomalies = process.shear_wave_to_density(dvs, zero_shallow_mantle=True, depths=depths)

latitudes = latitudes[:-1]
longitudes = np.roll(longitudes, shift=180)
density_anomalies = np.array(
    [np.roll(density_anomalies[i], shift=180, axis=1) for i in range(density_anomalies.shape[0])])[:, :-1, :]

# Convert to clm object for each depth.
clms = np.array([pysh.SHCoeffs.from_array(
    pysh.expand.SHExpandDH(griddh=density_anomalies[i], sampling=2, lmax_calc=l_max, norm=4, csphase=-1),
    normalization='ortho', csphase=-1) for i, _ in enumerate(depths)])


kernel_root = f'lyness/OUTPUT_{visc_src}/KERNELS/'
visc_name = visc_src

clm_DT_array = np.zeros((2, l_max + 1, l_max + 1))
density_anomalies_coeffs_array = np.array([clm.coeffs for clm in clms])

for degree in range(1, l_max + 1):
    kernel_path = kernel_root + kernel_type + read_in.leading_zeros(degree)

    radius_kernel, kernel = read_in.parse_file(kernel_path)
    radius_kernel *= 1.e3

    # Interpolate
    kernel = np.interp(depths, radius_kernel, kernel)

    integral = process.calculate_observable_at_degree(observable=kernel_type, l=degree, l_max=l_max,
                                            density_anomaly_sh_lm=density_anomalies_coeffs_array.real,
                                            kernel=kernel,
                                            radius_arr=depths)

    m_negative = integral[:l_max + 1]
    m_positive = integral[l_max + 1:]

    # Positive
    clm_DT_array[0][degree] = m_negative
    # Negative
    clm_DT_array[1][degree] = m_positive

viscosity = process.get_viscosity(f'lyness/VISC_INPUTS/{visc_src}.vis')
viscosity_radius = radius_kernel / 1.e3
fig, ax = plt.subplots(figsize=(6, 10))
ax.plot(viscosity, viscosity_radius, color='black')
ax.set_xscale('log')
ax.set_xlabel('$\eta/\eta_0$')
ax.set_ylabel('Radius (km)')
ax.set_ylim(viscosity_radius.min(), viscosity_radius.max())
# ax.set_xlim([0.5, 5e4])
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

plt.savefig(f'figures/S20RTS/{visc_src}/viscosity_profile.png', dpi=300)

# Plot DT
clm_DT = pysh.SHCoeffs.from_array(clm_DT_array, normalization='ortho')
clm_DT.csphase = -1

# Plot power spectrum
DT_power = process.plot_power_spectrum(clm_DT, unit='km2', save=True, fname=f'figures/S20RTS/{visc_src}/{kernel_type}/power_spectrum')


# plt.savefig('../../figures/S20RTS/DT_l2norm_.png', dpi=300)



# Plot dynamic topography
latitudes = np.linspace(-90, 90, 181)
longitudes = np.linspace(0, 360, 361)

longs, lats = np.meshgrid(longitudes, latitudes)

DT_grid = []
for lat in tqdm(lats):
    DT_grid.append(clm_DT.expand(grid='DH2', lat=lat, lon=longitudes))

DT_grid = np.array(DT_grid) / 1.e3

fig, ax = plt.subplots()
plt.imshow(DT_grid, cmap='RdBu_r', extent=[0, 360, -90, 90], origin='lower')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(extend='both', shrink=0.5, label=f'{kernel_type} (km)')
# plt.savefig('../../figures/S20RTS/DT_no_continents.png', dpi=300)


# Plot dynamic topography with continents
xr_DT = process.np_to_xr(DT_grid, lats=latitudes, lons=longitudes)

fig = pygmt.Figure()
fig.grdimage(xr_DT, cmap='jet', shading=True, projection='H12c', frame=True)
fig.coast(shorelines=True, resolution='c', area_thresh=1000, projection='H12c', region=[0, 360, -90, 90])
fig.colorbar(
    # Label the x-axis "Velocity" and the y-axis "m/s"
    frame=[f"x+l{kernel_type}", "y+l(km)"],

)
print('Saving DT...')
fig.savefig(f'figures/S20RTS/{visc_src}/{kernel_type}/map.png', dpi=300)
print('Done. \n')

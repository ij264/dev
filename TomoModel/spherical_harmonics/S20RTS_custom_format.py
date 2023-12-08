import numpy as np
import pyshtools as pysh
import matplotlib.pyplot as plt
import process
import read_in
import pygmt

# Set the maximum spherical harmonic degree.
l_max = 20

# Read in spherical harmonic coefficients
file_path = '../out/output.txt'
v_s_anomaly_sh_coeffs, radius, ref_density, ref_shear_velocity = read_in.read_in_S20RTS(file_path)


# Convert shear wave velocity coefficients to density coefficients.
density_anomaly_sh_coeffs = process.shear_wave_to_density(v_s_anomaly_sh_coeffs.real)

# Initialise a pyshtools SHCoeffs class for each radius.
density_coeffs_clm = np.array([pysh.SHCoeffs.from_zeros(lmax=l_max, kind='real',
                                                        normalization='ortho') for i, v in
                               enumerate(radius)])

# Use the Condon-Shortley phase convention.
for i, clm in enumerate(density_coeffs_clm):
    clm.csphase = 1
    clm = process.S20RTS_to_pyshtools(sh_coeffs=density_anomaly_sh_coeffs[i, :], clm=clm)

density_coeffs_array = np.array([clm.coeffs for clm in density_coeffs_clm])

kernel_root = '../data/kernels/'
visc_name = 'const_visc/'
kernel_type = 'surftopo'

clm_DT_array = np.zeros((2, l_max + 1, l_max + 1))

for degree in range(1, l_max + 1):
    kernel_path = kernel_root + visc_name + kernel_type + read_in.leading_zeros(degree)

    radius_kernel, kernel = read_in.parse_file(kernel_path)
    radius_kernel *= 1.e3

    # Interpolate
    kernel = np.interp(radius, radius_kernel, kernel)

    integral = process.calculate_observable(observable='surftopo', l=degree, l_max=l_max,
                                            density_anomaly_sh_lm=density_coeffs_array.real,
                                            kernel=kernel,
                                            radius_arr=radius)

    m_negative = integral[:l_max + 1]
    m_positive = integral[l_max + 1:]

    # Positive
    clm_DT_array[0][degree] = m_negative
    # Negative
    clm_DT_array[1][degree] = m_positive


# Plot DT
clm_DT = pysh.SHCoeffs.from_array(clm_DT_array, normalization='ortho')
clm_DT.csphase = 1
grid = clm_DT.expand(grid='DH') / 1.e3
fig, ax = grid.plot(colorbar='right',
                    cb_label='Dynamic Topography (km)',
                    cmap='viridis')
plt.show()
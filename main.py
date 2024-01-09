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
if len(sys.argv) != 1:
    print("Usage: python main.py <viscosity> <kernel_type>")
    sys.exit(1)

kernel_types = ['geoid', 'surftopo']

readable_kernels = {
    'surftopo': 'Dynamic Topography',
    'geoid': 'Geoid Anomaly'
}

# visc_src = sys.argv[1]
# kernel_type = sys.argv[1]
viscosities_location = 'lyness/VISC_INPUTS/step_changes'
visc_srces = os.listdir(viscosities_location)
visc_srces = [visc_src[:-4] for visc_src in visc_srces if not '_2' in visc_src]

# assert isinstance(visc_src, str), "Kernel source must be a string."
# assert kernel_type in valid_kernels, f"Kernel type must be one of {valid_kernels}."

pygmt.config(COLOR_BACKGROUND="0.6375/0.6375/255", COLOR_FOREGROUND="255/0.6375/0.6375")

# Maximum degrees
l_maxes = [3, 20]

# Set debugging toggle
debugging = False
location = 'figures'
if debugging:
    location = 'debugging'

# Controls whether the top 400 km is removed or not
excise_top_layer = True
excised_top_layer_location = 'excised'
if not excise_top_layer:
    excised_top_layer_location = 'whole_mantle'

# Loop over kernels
for kernel_type in kernel_types:
    match kernel_type:
        case 'surftopo':
            interval = 250
            annotation = 500
            series = [-2000, 2000, 10]
        case 'geoid':
            interval = 25
            annotation = 50
            series = [-100, 100, 10]

    # Loop over the maximum degree
    for l_max in l_maxes:
        # Loop over the viscosity profiles
        for visc_src in tqdm(visc_srces):
            # Create kernel if it does not exist.
            kernel_location = f'lyness/step_changes/OUTPUT_{visc_src}/'
            vsc_location = f'lyness/VISC_INPUTS/step_changes/{visc_src}'
            file_save_location = f'{location}/S20RTS/l_max = {l_max}/{excised_top_layer_location}/step_changes/{visc_src}/{kernel_type}/'

            if not os.path.exists(kernel_location):
                command = f'cd lyness && ./MAKE_KERNEL {kernel_location[7:]} {vsc_location[7:]}.vis && cd ..'

                # Run the command in the shell
                try:
                    result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            text=True)
                    print(result.stdout)
                    print("Kernels produced successfully")
                except subprocess.CalledProcessError as e:
                    print(f"Command failed with error: {e}")
                    print(e.stderr)

            if not os.path.exists(file_save_location):
                os.makedirs(file_save_location)

            # Read in data
            depths, latitudes, longitudes, dvs = read_in.read_in_nc("data/S20RTS_dvs.nc")

            # Convert to percentage.
            dvs /= 100
            # Convert to metres.
            depths *= 1.e3

            # Convert to shear wave velocity to density.
            density_anomalies = process.shear_wave_to_density(dvs, zero_shallow_mantle=excise_top_layer, depths=depths)

            # Note that longitude is from -180 to 180, so we need to shift it to 0 to 360.
            latitudes = latitudes[:-1]
            longitudes = np.roll(longitudes, shift=180)
            density_anomalies = np.array(
                [np.roll(density_anomalies[i], shift=180, axis=1) for i in range(density_anomalies.shape[0])])[:, :-1,
                                :]

            # Convert to clm object for each depth.
            clms = np.array([pysh.SHCoeffs.from_array(
                pysh.expand.SHExpandDH(griddh=density_anomalies[i], sampling=2, lmax_calc=l_max, norm=4, csphase=-1),
                normalization='ortho', csphase=-1) for i, _ in enumerate(depths)])

            kernel_root = f'{kernel_location}/KERNELS/'
            visc_name = visc_src

            clm_DT_array = np.zeros((2, l_max + 1, l_max + 1))
            density_anomalies_coeffs_array = np.array([clm.coeffs for clm in clms])

            for degree in range(1, l_max + 1):
                kernel_path = kernel_root + kernel_type + read_in.leading_zeros(degree)

                radius_kernel, kernel = read_in.parse_file(kernel_path)
                radius_kernel *= 1.e3

                # Interpolate
                # Generating indices for the original array
                original_indices = np.arange(len(kernel))

                # Generating indices for the new array
                new_indices = np.linspace(0, len(kernel) - 1, len(depths))

                # Interpolating the array to size 20
                interp_kernel = np.interp(new_indices, original_indices, kernel)

                integral = process.calculate_observable_at_degree(observable=kernel_type, l=degree, l_max=l_max,
                                                                  density_anomaly_sh_lm=density_anomalies_coeffs_array.real,
                                                                  kernel=interp_kernel,
                                                                  radius_arr=depths)

                m_negative = integral[:l_max + 1]
                m_positive = integral[l_max + 1:]

                # Positive
                clm_DT_array[0][degree] = m_negative
                # Negative
                clm_DT_array[1][degree] = m_positive

            viscosity = process.get_viscosity(f'{vsc_location}.vis')
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

            # plt.savefig(f'{location}/S20RTS/l_max = {l_max}/{excised_top_layer_location}/{visc_src}/viscosity_profile.png', dpi=300)

            # Plot DT
            clm_DT = pysh.SHCoeffs.from_array(clm_DT_array, normalization='ortho')
            clm_DT.csphase = -1

            # Plot power spectrum
            DT_power = process.plot_power_spectrum(clm_DT, unit='km2', save=True,
                                                   fname=f'{file_save_location}/power_spectrum')

            # Plot dynamic topography
            latitudes = np.linspace(-90, 90, 181)
            longitudes = np.linspace(0, 360, 361)

            longs, lats = np.meshgrid(longitudes, latitudes)

            DT_grid = []
            for lat in lats:
                DT_grid.append(clm_DT.expand(grid='DH2', lat=lat, lon=longitudes))

            DT_grid = np.array(DT_grid)

            # Plot dynamic topography with continents
            xr_DT = process.np_to_xr(DT_grid, lats=latitudes, lons=longitudes)

            fig = pygmt.Figure()

            pygmt.makecpt(
                cmap="polar",
                series=series,
                continuous=True
            )

            with fig.subplot(
                    nrows=1, ncols=2, figsize=("27.5c", "4c")):
                # Plot the original digital elevation model in the first panel
                with fig.set_panel(panel=0):
                    fig.grdimage(xr_DT,
                                 shading=None,
                                 projection='H12c',
                                 frame=f"t+t{readable_kernels[kernel_type]}",
                                 interpolation='b')
                    fig.grdcontour(xr_DT,
                                   interval=interval,
                                   annotation=annotation,
                                   limit=[-2000, 2000],
                                   projection='H12c')
                    fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H12c',
                              region=[0, 360, -90, 90])

                    fig.colorbar(frame=["x+lElevation",
                                        "y+lm"])
                # elevation model
                with fig.set_panel(panel=1):
                    fig.plot(
                        projection='X3cl/4c',
                        x=viscosity,
                        y=viscosity_radius,
                        frame=['yaf+lRadius (km)',
                               "xaf+leta (Pa s)"
                               ],
                        pen='1p',
                        region=[viscosity.min() / 2, viscosity.max() * 2, viscosity_radius.min(),
                                viscosity_radius.max()]
                    )
            print(f'Saving {kernel_type}...')
            fig.savefig(
                f'{file_save_location}/map.png',
                dpi=300)
            # print('Done. \n')

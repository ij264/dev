import process
import numpy as np
import pyshtools as pysh
import read_in
from tqdm import tqdm
import pygmt
import sys
import subprocess
import os

kernel_types = ['geoid', 'surftopo']

readable_kernels = {
    'surftopo': 'Dynamic Topography',
    'geoid': 'Geoid Anomaly'
}

viscosities_location = 'lyness/VISC_INPUTS/'
visc_srces = os.listdir(viscosities_location)
visc_srces = [visc_src[:-4] for visc_src in visc_srces if not '_2' in visc_src and not 'step_changes' in visc_src]

pygmt.config(COLOR_BACKGROUND="0.6375/0.6375/255", COLOR_FOREGROUND="255/0.6375/0.6375")

# Maximum degrees
l_maxes = [3, 20]

# Set debugging toggle
debugging = False
location = 'figures'
if debugging:
    location = 'debugging'

# Controls whether the top 400 km is removed or not
excise_top_layer = False
excised_top_layer_location = 'excised'
if not excise_top_layer:
    excised_top_layer_location = 'whole_mantle'

tomographic_model_root = 'coefficients/vs_coefficients'
tomographic_models = [model for model in os.listdir(tomographic_model_root)]

for visc_src in visc_srces:
    kernel_location = f'lyness/OUTPUT_{visc_src}/'
    vsc_location = f'lyness/VISC_INPUTS/{visc_src}'

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


def get_degree(str):
    return int(str[4:-11])


for tomographic_model in tomographic_models:

    files = os.listdir(f"{tomographic_model_root}/{tomographic_model}")
    files.sort(key=get_degree)
    data = np.load(f"{tomographic_model_root}/{tomographic_model}/{files[-1]}")
    clm_arrays = data['clm_arrays']
    depths = data['depths']
    latitudes = data['lats']
    longitudes = data['longs']

    cutoff_index_upper = np.argmax(depths >= 670)
    cutoff_index_lower = np.argmax(depths >= 100)

    upper_average_dvs, lower_average_dvs = process.average_shear_wave(clm_arrays,
                                                                      cutoff_index_lower=cutoff_index_lower,
                                                                      cutoff_index_upper=cutoff_index_upper)

    upper_average_dvs_arr = pysh.expand.MakeGrid2D(cilm=upper_average_dvs,
                                                   interval=1,
                                                   norm=4,
                                                   csphase=-1,
                                                   north=89,
                                                   south=-90,
                                                   west=0,
                                                   east=359
                                                   )

    lower_average_dvs_arr = pysh.expand.MakeGrid2D(cilm=lower_average_dvs,
                                                   interval=1,
                                                   norm=4,
                                                   csphase=-1,
                                                   north=89,
                                                   south=-90,
                                                   west=0,
                                                   east=359
                                                   )

    dvs_params = {
        'cmap': 'polar_dvs.cpt',
        'interval': 0.5,
        'annotation': 1,
    }

    xr_upper_dvs = process.np_to_xr(upper_average_dvs_arr, lats=latitudes, lons=longitudes)
    xr_lower_dvs = process.np_to_xr(lower_average_dvs_arr, lats=latitudes, lons=longitudes)

    surface_dvs = pysh.expand.MakeGrid2D(cilm=clm_arrays[0],
                                         interval=1,
                                         norm=4,
                                         csphase=-1,
                                         north=89,
                                         south=-90,
                                         west=0,
                                         east=359
                                         )

    surface_rho_clm = process.shear_wave_to_density(clm_arrays, depths, zero_shallow_mantle=False)[0]

    surface_rho = pysh.expand.MakeGrid2D(cilm=surface_rho_clm,
                                         interval=1,
                                         norm=4,
                                         csphase=-1,
                                         north=89,
                                         south=-90,
                                         west=0,
                                         east=359
                                         )

    # Plot dynamic topography with continents
    xr_surface_rho = process.np_to_xr(surface_rho, lats=latitudes, lons=longitudes)

    fig = pygmt.Figure()

    with fig.subplot(
            nrows=1,
            ncols=1,
            figsize=("27.5c", "16")
    ):
        with fig.set_panel(panel=0):
            fig.grdimage(xr_surface_rho,
                         cmap='polar_dvs.cpt',
                         shading=None,
                         projection='H10c',
                         frame=["a"],
                         interpolation='b')
            # fig.grdcontour(xr_surface_dvs,
            #                interval=1.,
            #                projection='H10c',
            #                frame=True)
            fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H10c',
                      region="d")

            fig.colorbar(position="+w12c/0.4c+h+o-1c/-1.8c",
                         cmap='polar_dvs.cpt',
                         frame=["y+l%"])
    fig.savefig(f'debugging/checking_conversion/rho.png')

    for l_max in l_maxes:
        # Read in tomographic_models
        if not os.path.exists(f"{tomographic_model_root}/{tomographic_model}/l = {l_max} coeffs.npz"):
            print(f"No l_max of {l_max} found. Moving onto the next one.")
            continue

        data = np.load(f"{tomographic_model_root}/{tomographic_model}/l = {l_max} coeffs.npz")
        clm_arrays = data['clm_arrays']
        # Convert to percentage.
        clm_arrays /= 100
        # Convert to metres.
        depths *= 1.e3
        # Convert to shear wave velocity to density.
        density_anomalies_clm = process.shear_wave_to_density(clm_arrays, zero_shallow_mantle=excise_top_layer,
                                                              depths=depths)

        for kernel_type in kernel_types:
            print(f'\n\nTomographic model: {tomographic_model} \nl max = {l_max} \nkernel type: {kernel_type}')
            match kernel_type:
                case 'surftopo':
                    interval = 0.2
                    annotation = 0.5
                    series = [-2000, 2000, 10]
                    scale = 1.e3
                    cpt = 'polar_DT.cpt'
                    limits = [-2., 2.]
                    unit = 'km'
                case 'geoid':
                    interval = 10
                    annotation = 50
                    series = [-100, 100, 10]
                    scale = 1
                    cpt = 'polar_geoid.cpt'
                    limits = [-200, 200]
                    unit = 'm'

            for visc_src in tqdm(visc_srces):
                # Create kernel if it does not exist.
                kernel_location = f'lyness/OUTPUT_{visc_src}'
                vsc_location = f'lyness/VISC_INPUTS/{visc_src}'
                file_save_location = f'{location}/{tomographic_model}/l_max = {l_max}/{excised_top_layer_location}/{visc_src}/{kernel_type}'

                if not os.path.exists(file_save_location):
                    os.makedirs(file_save_location)
                kernel_root = f'{kernel_location}/KERNELS/'
                visc_name = visc_src

                clm_observable_array = np.zeros((2, l_max + 1, l_max + 1))

                # Integrate over the depth
                for degree in range(1, l_max + 1):
                    kernel_path = kernel_root + kernel_type + read_in.leading_zeros(degree)

                    radius_kernel, kernel = read_in.parse_file(kernel_path)
                    radius_kernel *= 1.e3

                    # Interpolate the kernel and radius to fit the tomographic data
                    interp_kernel = process.interp(kernel, len(depths))

                    # Calculate the integral
                    integral = process.calculate_observable_at_degree(observable=kernel_type, l=degree, l_max=l_max,
                                                                      density_anomaly_sh_lm=density_anomalies_clm,
                                                                      kernel=interp_kernel,
                                                                      radius_arr=depths)

                    m_negative = integral[:l_max + 1]
                    m_positive = integral[l_max + 1:]

                    # Positive
                    clm_observable_array[0][degree] = m_positive
                    # Negative
                    clm_observable_array[1][degree] = m_negative

                # Plot viscosity
                viscosity = process.get_viscosity(f'{vsc_location}.vis')
                viscosity_radius = radius_kernel / 1.e3

                # Save coefficients
                if not os.path.exists(f'{tomographic_model_root}/{kernel_type}/coefficients/{visc_name}'):
                    os.makedirs(f'{tomographic_model_root}/{kernel_type}/coefficients/{visc_name}')
                pysh.shio.shwrite(
                    filename=f'{tomographic_model_root}/{kernel_type}/coefficients/{visc_name}/{kernel_type} lmax = {l_max}.coef',
                    coeffs=clm_observable_array)

                # Plot DT
                clm_observable = pysh.SHCoeffs.from_array(clm_observable_array, normalization='ortho')
                clm_observable.csphase = -1

                # Plot power spectrum

                observable_power = process.plot_power_spectrum(clm_observable, unit='km2', save=True,
                                                               fname=f'{file_save_location}/power_spectrum')

                # Plot dynamic topography
                observable_grid = pysh.expand.MakeGrid2D(cilm=clm_observable_array,
                                                         interval=1,
                                                         norm=4,
                                                         csphase=-1,
                                                         north=89,
                                                         south=-90,
                                                         west=0,
                                                         east=359
                                                         )
                observable_grid /= scale

                # DT_grid = np.roll(DT_grid, shift=180, axis=1) / scale
                # DT_grid /= 1.e3

                # Plot dynamic topography with continents
                xr_observable = process.np_to_xr(observable_grid / scale, lats=latitudes, lons=longitudes)

                # Plot dynamic topography with continents

                fig = pygmt.Figure()

                with fig.subplot(
                        nrows=2,
                        ncols=2,
                        figsize=("27.5c", "16")
                ):
                    # Plot the original digital elevation model in the first panel
                    with fig.set_panel(panel=0):
                        fig.grdimage(xr_observable,
                                     cmap=cpt,
                                     shading=None,
                                     projection='H10c',
                                     frame=["a", f"t+t{readable_kernels[kernel_type]}"],
                                     interpolation='b')
                        fig.grdcontour(xr_observable,
                                       interval=interval,
                                       # annotation=annotation,
                                       projection='H10c',
                                       frame=True)
                        fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H10c',
                                  region="d")

                        fig.colorbar(position="+w12c/0.4c+h+o-1c/-1.8c",
                                     cmap=cpt,
                                     frame=[f"y+l{unit}"],
                                     panel=0)

                    with fig.set_panel(panel=1):
                        if np.any(viscosity <= 0):
                            proj = 'X4c/6c'
                        else:
                            proj = 'X4cl/6c'
                        fig.plot(
                            projection=proj,
                            x=viscosity,
                            y=viscosity_radius,
                            frame=['yaf+lRadius (km)',
                                   "xaf+leta (Pa s)"
                                   ],
                            pen='1p',
                            region=[viscosity.min() / 2, viscosity.max() * 2, viscosity_radius.min(),
                                    viscosity_radius.max()],
                            panel=1
                        )
                    with fig.set_panel(panel=2):
                        fig.grdimage(xr_upper_dvs,
                                     cmap=dvs_params['cmap'],
                                     shading=None,
                                     projection='H8c',
                                     frame=["a", "t+tUpper mantle"],
                                     interpolation='b')
                        fig.grdcontour(xr_upper_dvs,
                                       interval=dvs_params['interval'],
                                       projection='H8c',
                                       frame=True)
                        fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H8c',
                                  region="d")

                        fig.colorbar(position="+w10c/0.25c+h+o-1c/-1.5c",
                                     cmap=dvs_params['cmap'],
                                     frame=["x+ldv/v", "y+l%"],
                                     panel=2)

                    with fig.set_panel(panel=3):
                        fig.grdimage(xr_lower_dvs,
                                     cmap=dvs_params['cmap'],
                                     shading=None,
                                     projection='H8c',
                                     frame=["a", "t+tLower mantle"],
                                     interpolation='b')
                        fig.grdcontour(xr_lower_dvs,
                                       interval=interval,
                                       annotation=annotation,
                                       projection='H8c',
                                       frame=True)
                        fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H8c',
                                  region="d")

                        fig.colorbar(position="+w10c/0.25c+h+o-1c/-1.5c",
                                     cmap=dvs_params['cmap'],
                                     frame=["x+ldv/v", "y+l%"],
                                     panel=3
                                     )

                print(f'Saving {kernel_type}...')
                fig.savefig(
                    f'{file_save_location}/map.png',
                    dpi=300)
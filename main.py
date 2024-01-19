# TODO
# Figure out why it breaks if there are multiple l_maxes. Seems to break down at the integral for degree = 1 of the second l_max. Integral is FAR too large.


import process
import numpy as np
import pyshtools as pysh
import read_in
from tqdm import tqdm
import pygmt
import sys
import subprocess
import os
import matplotlib.pyplot as plt

# Retrieve command-line arguments
if len(sys.argv) != 1:
    print("Usage: python main.py <viscosity> <kernel_type>")
    sys.exit(1)

# kernel_types = ['geoid', 'surftopo']
kernel_types = ['geoid', 'surftopo']
readable_kernels = {
    'surftopo': 'Dynamic Topography',
    'geoid': 'Geoid Anomaly'
}

# visc_src = sys.argv[1]
# kernel_type = sys.argv[1]
viscosities_location = 'lyness/VISC_INPUTS/'
visc_srces = os.listdir(viscosities_location)
# visc_srces = [visc_src[:-4] for visc_src in visc_srces if not '_2' in visc_src]
visc_srces = [visc_src[:-4] for visc_src in visc_srces if not '_2' in visc_src and not 'step_changes' in visc_src]
# assert isinstance(visc_src, str), "Kernel source must be a string."
# assert kernel_type in valid_kernels, f"Kernel type must be one of {valid_kernels}."

pygmt.config(COLOR_BACKGROUND="0.6375/0.6375/255", COLOR_FOREGROUND="255/0.6375/0.6375")

# Maximum degrees
l_maxes = [3, 20]
# l_max = 40
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
# tomographic_models = [model for model in os.listdir(tomographic_model_root) if 'S40' in model]

for visc_src in visc_srces:
    # Create kernel if it does not exist.
    kernel_location = f'lyness/OUTPUT_{visc_src}/'
    vsc_location = f'lyness/VISC_INPUTS/{visc_src}'
    # file_save_location = f'{location}/{tomographic_model[:-7]}/l_max = {l_max}/{excised_top_layer_location}/{visc_src}/{kernel_type}/'

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
#
# for tomographic_model in tqdm(tomographic_models):
#     # Read in tomographic_models
#     files = os.listdir(f"{tomographic_model_root}/{tomographic_model}")
#     files.sort(key=get_degree)
#     data = np.load(f"{tomographic_model_root}/{tomographic_model}/{files[-1]}")
#     clm_arrays = data['clm_arrays']
#     depths = data['depths']
#     latitudes = data['lats']
#     longitudes = data['longs']
#
#     surface_dvs = pysh.expand.MakeGrid2D(cilm=clm_arrays[0],
#                                          interval=1,
#                                          norm=4,
#                                          csphase=-1,
#                                          north=89,
#                                          south=-90,
#                                          west=0,
#                                          east=359
#                                          )
#     xr_surface_dvs = process.np_to_xr(surface_dvs, lats=latitudes, lons=longitudes)
#
#     # surface_dvs = np.roll(surface_dvs, shift=180, axis=1)
#     # DT_grid /= 1.e3
#
#     # Plot dynamic topography with continents
#     xr_surface_dvs = process.np_to_xr(surface_dvs, lats=latitudes, lons=longitudes)
#
#     fig = pygmt.Figure()
#
#     with fig.subplot(
#             nrows=1,
#             ncols=1,
#             figsize=("27.5c", "16")
#     ):
#         # Plot the original digital elevation model in the first panel
#         with fig.set_panel(panel=0):
#             fig.grdimage(xr_surface_dvs,
#                          cmap='polar_dvs.cpt',
#                          shading=None,
#                          projection='H10c',
#                          frame=["a"],
#                          interpolation='b')
#             fig.grdcontour(xr_surface_dvs,
#                            interval=1.,
#                            projection='H10c',
#                            frame=True)
#             fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H10c',
#                       region="d")
#
#             fig.colorbar(position="+w12c/0.4c+h+o-1c/-1.8c",
#                          cmap='polar_dvs.cpt',
#                          frame=["y+l%"])
#     print(tomographic_model)
#     fig.savefig(f'debugging/testing_tomography_models/{tomographic_model}.png')


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

    # upper_average_dvs_arr = np.roll(upper_average_dvs_arr, shift=180, axis=1)
    # lower_average_dvs_arr = np.roll(lower_average_dvs_arr, shift=180, axis=1)

    xr_upper_dvs = process.np_to_xr(upper_average_dvs_arr, lats=latitudes, lons=longitudes)
    xr_lower_dvs = process.np_to_xr(lower_average_dvs_arr, lats=latitudes, lons=longitudes)

    # surface_dvs = pysh.expand.MakeGrid2D(cilm=clm_arrays[0],
    #                                      interval=1,
    #                                      norm=4,
    #                                      csphase=-1,
    #                                      north=89,
    #                                      south=-90,
    #                                      west=0,
    #                                      east=359
    #                                      )
    # xr_surface_dvs = process.np_to_xr(surface_dvs, lats=latitudes, lons=longitudes)
    #
    # # surface_dvs = np.roll(surface_dvs, shift=180, axis=1)
    # # DT_grid /= 1.e3
    #
    # # Plot dynamic topography with continents
    # xr_surface_dvs = process.np_to_xr(surface_dvs, lats=latitudes, lons=longitudes)
    #
    # fig = pygmt.Figure()
    #
    # with fig.subplot(
    #         nrows=1,
    #         ncols=1,
    #         figsize=("27.5c", "16")
    # ):
    #     # Plot the original digital elevation model in the first panel
    #     with fig.set_panel(panel=0):
    #         fig.grdimage(xr_surface_dvs,
    #                      cmap='polar_dvs.cpt',
    #                      shading=None,
    #                      projection='H10c',
    #                      frame=["a"],
    #                      interpolation='b')
    #         # fig.grdcontour(xr_surface_dvs,
    #         #                interval=1.,
    #         #                projection='H10c',
    #         #                frame=True)
    #         fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H10c',
    #                   region="d")
    #
    #         fig.colorbar(position="+w12c/0.4c+h+o-1c/-1.8c",
    #                      cmap='polar_dvs.cpt',
    #                      frame=["y+l%"])
    # fig.savefig(f'debugging/checking_conversion/v_s.png')

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

    # surface_dvs = np.roll(surface_dvs, shift=180, axis=1)
    # DT_grid /= 1.e3
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
        # Plot the original digital elevation model in the first panel
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

                clm_DT_array = np.zeros((2, l_max + 1, l_max + 1))

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
                    clm_DT_array[0][degree] = m_positive
                    # Negative
                    clm_DT_array[1][degree] = m_negative

                # Plot viscosity
                viscosity = process.get_viscosity(f'{vsc_location}.vis')
                viscosity_radius = radius_kernel / 1.e3

                # Save coefficients
                if not os.path.exists(f'{tomographic_model_root}/{kernel_type}/coefficients/{visc_name}'):
                    os.makedirs(f'{tomographic_model_root}/{kernel_type}/coefficients/{visc_name}')
                pysh.shio.shwrite(
                    filename=f'{tomographic_model_root}/{kernel_type}/coefficients/{visc_name}/{kernel_type} lmax = {l_max}.coef',
                    coeffs=clm_DT_array)

                # Plot DT
                clm_DT = pysh.SHCoeffs.from_array(clm_DT_array, normalization='ortho')
                clm_DT.csphase = -1

                # Plot power spectrum

                DT_power = process.plot_power_spectrum(clm_DT, unit='km2', save=True,
                                                       fname=f'{file_save_location}/power_spectrum')

                # Plot dynamic topography
                DT_grid = pysh.expand.MakeGrid2D(cilm=clm_DT_array,
                                                 interval=1,
                                                 norm=4,
                                                 csphase=-1,
                                                 north=89,
                                                 south=-90,
                                                 west=0,
                                                 east=359
                                                 )
                DT_grid /= scale

                # DT_grid = np.roll(DT_grid, shift=180, axis=1) / scale
                # DT_grid /= 1.e3

                # Plot dynamic topography with continents
                xr_DT = process.np_to_xr(DT_grid / scale, lats=latitudes, lons=longitudes)

                # Plot dynamic topography with continents

                fig = pygmt.Figure()

                with fig.subplot(
                        nrows=2,
                        ncols=2,
                        figsize=("27.5c", "16")
                ):
                    # Plot the original digital elevation model in the first panel
                    with fig.set_panel(panel=0):
                        fig.grdimage(xr_DT,
                                     cmap=cpt,
                                     shading=None,
                                     projection='H10c',
                                     frame=["a", f"t+t{readable_kernels[kernel_type]}"],
                                     interpolation='b')
                        fig.grdcontour(xr_DT,
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

                # with fig.subplot(
                #         nrows=1, ncols=2, figsize=("27.5c", "4c")):
                #     # Plot the original digital elevation model in the first panel
                #     with fig.set_panel(panel=0):
                #         fig.grdimage(xr_DT,
                #                      cmap='coolwarm_DT.cpt',
                #                      shading=None,
                #                      projection='H12c',
                #                      frame=["a", f"t+t{readable_kernels[kernel_type]}"],
                #                      interpolation='b')
                #         fig.grdcontour(xr_DT,
                #                        interval=interval,
                #                        annotation=annotation,
                #                        limit=[-2.000, 2.000],
                #                        projection='H12c',
                #                        frame=True)
                #         fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H12c',
                #                   region="d")
                #
                #         fig.colorbar(cmap='coolwarm_DT.cpt',
                #                      frame=["x+lElevation",
                #                             "y+lm"])
                #     # elevation model
                #     with fig.set_panel(panel=1):
                #         fig.plot(
                #             projection='X3cl/4c',
                #             x=viscosity,
                #             y=viscosity_radius,
                #             frame=['yaf+lRadius (km)',
                #                    "xaf+leta (Pa s)"
                #                    ],
                #             pen='1p',
                #             region=[viscosity.min() / 2, viscosity.max() * 2, viscosity_radius.min(),
                #                     viscosity_radius.max()]
                #         )
                # print(f'Saving {kernel_type}...')
                # fig.savefig(
                #     f'{file_save_location}/map.png',
                #     dpi=300)
                # print('Done. \n')

        # # Note that longitude is from -180 to 180, so we need to shift it to 0 to 360.
        # # latitudes = latitudes[:-1]
        #
        # clm_arrays = []
        # LAT, LON = np.meshgrid(latitudes, longitudes, indexing='ij')
        # for i, depth in tqdm(enumerate(depths)):
        #     xx, yy, zz = LON.flatten(), LAT.flatten(), dvs[i]
        #     clm_arrays.append(
        #         [
        #             pysh.expand.SHExpandLSQ(d=zz, lat=xx, lon=yy, lmax=l_max)
        #         ]
        #     )
        # clm_arrays = np.array(clm_arrays)
        # np.savez_compressed(file=f'coefficients/{tomographic_model}/lmax = {l_max} coeffs.npz', lat=latitudes,
        #                     lon=longitudes, dvs=dvs)

#
# for tomographic_model in tomographic_models:
#     # Loop over kernels
#     for kernel_type in kernel_types:
#         match kernel_type:
#             case 'surftopo':
#                 interval = 250
#                 annotation = 500
#                 series = [-2000, 2000, 10]
#             case 'geoid':
#                 interval = 25
#                 annotation = 50
#                 series = [-100, 100, 10]
#
#         # Loop over the maximum degree
#         for l_max in l_maxes:
#             # Loop over the viscosity profiles
#             for visc_src in tqdm(visc_srces):
#                 # Create kernel if it does not exist.
#                 kernel_location = f'lyness/OUTPUT_{visc_src}/'
#                 vsc_location = f'lyness/VISC_INPUTS/{visc_src}'
#                 file_save_location = f'{location}/{tomographic_model[:-7]}/l_max = {l_max}/{excised_top_layer_location}/{visc_src}/{kernel_type}/'
#
#                 if not os.path.exists(kernel_location):
#                     command = f'cd lyness && ./MAKE_KERNEL {kernel_location[7:]} {vsc_location[7:]}.vis && cd ..'
#
#                     # Run the command in the shell
#                     try:
#                         result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE,
#                                                 stderr=subprocess.PIPE,
#                                                 text=True)
#                         print(result.stdout)
#                         print("Kernels produced successfully")
#                     except subprocess.CalledProcessError as e:
#                         print(f"Command failed with error: {e}")
#                         print(e.stderr)
#
#                 if not os.path.exists(file_save_location):
#                     os.makedirs(file_save_location)
#
#                 # Read in tomographic_models
#                 depths, latitudes, longitudes, dvs = read_in.read_in_nc(f"seismic-tomography-models/{tomographic_model}")
#
#                 # Convert to percentage.
#                 dvs /= 100
#                 # Convert to metres.
#                 depths *= 1.e3
#
#                 # Convert to shear wave velocity to density.
#                 density_anomalies = process.shear_wave_to_density(dvs, zero_shallow_mantle=excise_top_layer, depths=depths)
#
#                 # Note that longitude is from -180 to 180, so we need to shift it to 0 to 360.
#                 # latitudes = latitudes[:-1]
#
#                 clm_arrays = []
#                 LAT, LON = np.meshgrid(latitudes, longitudes, indexing='ij')
#                 for i, depth in tqdm(enumerate(depths)):
#                     xx, yy, zz = LON.flatten(), LAT.flatten(), density_anomalies[i]
#                     clm_arrays.append(
#                         [
#                         pysh.expand.SHExpandLSQ(d=zz, lat=xx, lon=yy, lmax=l_max)
#                         ]
#                     )
#                 clm_arrays = np.array(clm_arrays)
#
#
#
#
#
#
#
#                 # Convert to clm object for each depth.
#                 # clm_arrays = np.array([pysh.expand.SHExpandDH(griddh=density_anomalies[i], sampling=2, lmax_calc=l_max, norm=4, csphase=-1) for i, _ in enumerate(depths)])
#                 clms = [pysh.SHCoeffs.from_array(clm_array, normalization='ortho', csphase=-1) for clm_array in clm_arrays]
#
#                 kernel_root = f'{kernel_location}/KERNELS/'
#                 visc_name = visc_src
#
#                 clm_DT_array = np.zeros((2, l_max+1, l_max+1))
#
#                 # Integrate over the depth
#                 for degree in range(1, l_max + 1):
#                     kernel_path = kernel_root + kernel_type + read_in.leading_zeros(degree)
#
#                     radius_kernel, kernel = read_in.parse_file(kernel_path)
#                     radius_kernel *= 1.e3
#
#                     # Interpolate the kernel and radius to fit the tomographic data
#                     original_indices = np.arange(len(kernel))
#                     new_indices = np.linspace(0, len(kernel) - 1, len(depths))
#                     interp_kernel = np.interp(new_indices, original_indices, kernel)
#
#                     # Calculate the integral
#                     integral = process.calculate_observable_at_degree(observable=kernel_type, l=degree, l_max=l_max,
#                                                                       density_anomaly_sh_lm=clm_arrays,
#                                                                       kernel=interp_kernel,
#                                                                       radius_arr=depths)
#
#                     m_negative = integral[:l_max + 1]
#                     m_positive = integral[l_max + 1:]
#
#                     # Positive
#                     clm_DT_array[0][degree] = m_positive
#                     # Negative
#                     clm_DT_array[1][degree] = m_negative
#
#                 # Plot viscosity
#                 viscosity = process.get_viscosity(f'{vsc_location}.vis')
#                 viscosity_radius = radius_kernel / 1.e3
#
#                 # Save coefficients
#                 if not os.path.exists(f'{tomographic_model_root}/coefficients/{visc_name}'):
#                     os.makedirs(f'{tomographic_model_root}/coefficients/{visc_name}')
#                 pysh.shio.shwrite(filename=f'{tomographic_model_root}/coefficients/{visc_name}/{kernel_type} lmax = {l_max}.coef',
#                                   coeffs=clm_DT_array)
#
#                 # Plot DT
#                 clm_DT = pysh.SHCoeffs.from_array(clm_DT_array, normalization='ortho')
#                 clm_DT.csphase = -1
#
#                 # Plot power spectrum
#                 DT_power = process.plot_power_spectrum(clm_DT, unit='km2', save=True,
#                                                        fname=f'{file_save_location}/power_spectrum')
#
#                 # Plot dynamic topography
#                 latitudes = np.linspace(90, -90, 181)
#                 longitudes = np.linspace(0, 360, 361)
#
#                 DT_grid = pysh.expand.MakeGrid2D(cilm=clm_DT_array,
#                                                  interval=1,
#                                                  norm=4,
#                                                  csphase=-1)
#
#                 DT_grid = np.roll(DT_grid, shift=180, axis=1) / 1.e3
#
#                 # Plot dynamic topography with continents
#                 xr_DT = process.np_to_xr(DT_grid, lats=latitudes, lons=longitudes)
#
#                 fig = pygmt.Figure()
#
#                 with fig.subplot(
#                         nrows=1, ncols=2, figsize=("27.5c", "4c")):
#                     # Plot the original digital elevation model in the first panel
#                     with fig.set_panel(panel=0):
#                         fig.grdimage(xr_DT,
#                                      cmap='coolwarm_DT.cpt',
#                                      shading=None,
#                                      projection='H12c',
#                                      frame=["a", f"t+t{readable_kernels[kernel_type]}"],
#                                      interpolation='b')
#                         fig.grdcontour(xr_DT,
#                                        interval=interval,
#                                        annotation=annotation,
#                                        limit=[-2.000, 2.000],
#                                        projection='H12c',
#                                        frame=True)
#                         fig.coast(shorelines='1/0.75p,black', resolution='c', area_thresh=1000, projection='H12c',
#                                   region="d")
#
#                         fig.colorbar(cmap='coolwarm_DT.cpt',
#                                      frame=["x+lElevation",
#                                             "y+lm"])
#                     # elevation model
#                     with fig.set_panel(panel=1):
#                         fig.plot(
#                             projection='X3cl/4c',
#                             x=viscosity,
#                             y=viscosity_radius,
#                             frame=['yaf+lRadius (km)',
#                                    "xaf+leta (Pa s)"
#                                    ],
#                             pen='1p',
#                             region=[viscosity.min() / 2, viscosity.max() * 2, viscosity_radius.min(),
#                                     viscosity_radius.max()]
#                         )
#                 print(f'Saving {kernel_type}...')
#                 fig.savefig(
#                     f'{file_save_location}/map.png',
#                     dpi=300)
#                 # print('Done. \n')

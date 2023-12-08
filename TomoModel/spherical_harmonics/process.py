import pyshtools as pysh
import numpy as np
import scipy
from validation import assert_shape
from typing import Tuple
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def S20RTS_to_pyshtools(clm: pysh.SHCoeffs, sh_coeffs: np.ndarray) -> pysh.SHCoeffs:
    """
    Takes in a 1D NumPy array of S20RTS spherical harmonic coefficients and returns a 3D NumPy array of pyshtools spherical harmonic coefficents.

    Parameters
    ----------
    clm : pysh.SHCoeffs(cphase=1, normalization='ortho')
        Non-initialised spherical harmonic coefficient class.
    sh_coeffs : np.ndarray(dtype=float, ndim=1)
        S20RTS spherical harmonic coefficients.

    Returns
    -------
    clm : pysh.SHCoeffs(cphase=1, normalization='ortho')
        Initialised spherical harmonic coefficient class.
    """

    '''
    For positive orders, the coefficients are stored in the following order:

    c_{0, 0}           0,           0,       ...,       0
    c_{1, 0}        c{1, 1},        0,       ...,       0
    c_{2, 0}        c{2, 1},     c{2, 2},    ...,       0
    c_{3, 0}        c{3, 1},     c{3, 2},    ...,       0
        .              .            .         .        .        
        .              .            .         .        .         
        .              .            .         .        .         
    c_{l_max, 0}, c{l_max, 1}, c{l_max, 2}, ..., c{l_max, l_max}


    For negative orders, the coefficients are stored in the following order:

    0       0,            0,        ...,          0
    0   c_{1, -1},        0,        ...,          0
    0   c_{2, -1},    c_{2, -2},    ...,          0
    0   c_{3, -1},    c_{3, -2},    ...,          0
    .      .              .           .           .        
    .      .              .           .           .         
    .      .              .           .           .         
    0 c_{l_max, -1}, c_{l_max, -2}, ..., c_{l_max, -l_max}

    '''

    assert_shape(sh_coeffs, ((clm.lmax + 1) ** 2,))
    # Positive orders
    for l in range(clm.lmax + 1):
        index = l * (l + 1)
        for order in range(l + 1):
            clm.coeffs[0, l, order] = sh_coeffs[index + order]

    # Negative orders
    for l in range(1, clm.lmax + 1):
        index = l * (l + 1)
        for order in range(1, l + 1):
            clm.coeffs[1, l, order] = sh_coeffs[index - order]

    return clm


def read_in_nc(nc_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    This function takes in a netCDF file and returns the data.

    Parameters
    ----------
    nc_file : str
        The netCDF file to read in.

    Returns
    -------
    depths : np.ndarray(dtype=float, ndim=1)
        The depths that the data is sampled at. Units: km.
    latitudes : np.ndarray(dtype=float, ndim=1)
        The latitudes that the data is sampled at. Units: degrees.
    longitudes : np.ndarray(dtype=float, ndim=1)
        The longitudes that the data is sampled at. Units: degrees.
    v : np.ndarray(dtype=float, ndim=3)
        v is a 3D array of shape (depths, latitudes, longitudes). It can represent several variables.
    """
    import netCDF4 as nc

    data = nc.Dataset(nc_file)

    depths = data.variables['depth'][:]
    latitudes = data.variables['latitude'][:]
    longitudes = data.variables['longitude'][:]
    v = data.variables['v'][:]

    return depths, latitudes, longitudes, v


def shear_wave_to_density(dlnv_s: np.ndarray, depths, zero_shallow_mantle: bool = True) -> np.ndarray:
    """
    This function takes in shear wave velocity anomaly coefficients and returns density anomaly
    coefficients.

    Note: d ln rho / d ln v_s ~ 0.1

    Parameters
    ----------
    shear_wave_coeffs : np.ndarray(dtype=float, ndim=2)
        Array of shear wave velocity anomaly coefficients.
    depths : np.ndarray(dtype=float, ndim=1)
        Array of depths.
    zero_shallow_mantle : bool (optional)
        Zeros the top 400 km of density. Default: True.

    Returns
    -------
    delta_density_sh_lm : np.ndarray(dtype=float, ndim=2)
        Array of density anomaly coefficients using a scaling factor.
    """

    mean_density = 4500
    drho = 0.1 * dlnv_s * mean_density

    if zero_shallow_mantle:
        shallow_mantle_index = np.argmax(depths >= 400e3)
        drho[:shallow_mantle_index, :, :] = 0
    return drho


def calculate_observable_at_degree(observable: str, l: int, l_max: int, density_anomaly_sh_lm: np.ndarray,
                         radius_arr: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """
    This function calculates the observable for a given degree.

    The set of observables is:
        {
            'surftopo': Dynamic Topography,
            'geoid': Geoid,
            'gravity': Gravity Anomaly
            'cmbtopo': CMB Topography
        }

    Parameters
    ----------
    observable : str
        The observable to calculate.
    l : int
        The degree to calculate the observable at.
    l_max : int
        The maximum degree of the model.
    density_anomaly_sh_lm : np.ndarray(dtype=float, ndim=4)
        Array of density anomaly spherical harmonic coefficients.
    radius_arr : np.ndarray(dtype=float, ndim=1)
        Array of radii values.
    kernel : np.ndarray(dtype=float, ndim=1)
        Array of kernel values.

    Returns
    -------
    integral : np.ndarray(dtype=float, ndim=1)
        Array of integrated density anomalies x kernel.
    """
    assert observable in ['surftopo', 'geoid', 'gravity', 'cmbtopo'], f"Observable {observable} not recognised."

    rho_uppermost_mantle = 3380
    rho_water = 1030
    rho_lowermost_mantle = 5570
    rho_uppermost_outer_crust = 9900
    G = 6.67408e-11
    g = 9.81
    R = 6371e3

    match observable:
        case 'surftopo':
            init_factor = 1 / (rho_uppermost_mantle - rho_water)
        case 'geoid':
            init_factor = 4 * np.pi * G * R / ((2 * l + 1) * g)
        case 'gravity':
            init_factor = 1
        case 'cmbtopo':
            init_factor = -1 / (rho_lowermost_mantle - rho_uppermost_outer_crust)
        case _:
            raise ValueError(f"Observable {observable} not recognised.")

    assert_shape(density_anomaly_sh_lm, (len(radius_arr), 2, l_max + 1, l_max + 1))
    density_anomaly_sh_l = density_anomaly_sh_lm[:, :, l, :]

    # integrand is (radial_points x (2*l_max + 1))
    integrand = np.concatenate(
        (density_anomaly_sh_l[:, 1, :] * kernel.reshape(-1, 1),
         density_anomaly_sh_l[:, 0, :] * kernel.reshape(-1, 1)),
        axis=1
    )

    # integrate along the radius for each order
    integral = init_factor * scipy.integrate.simpson(integrand, radius_arr, axis=0)
    assert_shape(integral, (2 * (l_max + 1),))
    return integral


def np_to_xr(np_arr: np.ndarray, lats: np.ndarray, lons: np.ndarray) -> xr.DataArray:
    """
    Convert a 2D NumPy array to an xarray DataArray for use in PyGMT.

    Parameters
    ----------
    np_arr : np.ndarray
        2D NumPy array of values.
    lats : np.ndarray
        1D NumPy array of latitudes.
    lons : np.ndarray
        1D NumPy array of longitudes.

    Returns
    -------
    xr_arr : xr.DataArray
        xarray DataArray of values.
    """
    xr_arr = xr.DataArray(
        data=np_arr,
        coords=dict(
            y=lats,
            x=lons,
        ),
        dims=("y", "x"),
    )

    return xr_arr


def plot_power_spectrum(clm: pysh.SHCoeffs, unit='km2', save=False, **kwargs) -> np.ndarray:
    """
    This function takes in a pyshtools spherical harmonic coefficient class and returns the power spectrum.

    Parameters
    ----------
    clm : pysh.SHCoeffs(cphase=-1, normalization='ortho')
        Initialised spherical harmonic coefficient class.
    unit : str (optional)
        Values: ['km2, m2']
        The unit of the power spectrum. Default: km^2.
    save : bool (optional)
        Whether to save the plot. Default: False.
    **kwargs : dict
        Keyword arguments to pass to plt.subplots().

    Returns
    -------
    power_spectrum : np.ndarray(dtype=float, ndim=1)
        Array of power spectrum values.
    """
    _, ax = clm.plot_spectrum(convention='l2norm',
                              unit='per_lm',
                              show=False)

    power = ax.get_lines()[0].get_ydata()

    if unit:
        assert unit in ['km2', 'm2'], f"Unit {unit} not recognised."
        if unit == 'km2':
            power /= 1.e6

    fig, ax = plt.subplots()
    ax.plot(np.arange(1, clm.lmax + 1), power[1:], linestyle='--', marker='o', color='k')
    ax.set_xlabel('Degree, $l$')
    ax.set_xlim(1, clm.lmax)
    ax.set_ylim(1.e-2, 1.e1)
    ax.set_ylabel('Power [km$^2$]')
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    if save:
        print('Saving power spectrum plot...')
        plt.savefig(f"{kwargs['fname']}.png", dpi=300)
        print('Done.')


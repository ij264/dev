import matplotlib.pyplot as plt
import numpy as np
import pyshtools as pysh
from matplotlib.widgets import Slider


def S20RTS_to_pyshtools(clm, sh_coeffs):
    """
    Takes in a 1D NumPy array of S20RTS spherical harmonic coefficients and returns a 3D NumPy array of pyshtools spherical harmonic coefficents.


    :param pysh.SHCoeffs clm: Non-initialised spherical harmonic coefficient class.
    :param np.ndarray sh_coeffs: S20RTS spherical harmonic coefficients.
    :return pysh.SHCoeffs clm: Initialised spherical harmonic coefficient class.

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

    # Positive orders
    for l in range(l_max + 1):
        index = l * (l + 1)
        for order in range(l + 1):
            clm.coeffs[0, l, order] = sh_coeffs[index + order]

    # Negative orders
    for l in range(1, l_max + 1):
        index = l * (l + 1)
        for order in range(1, l + 1):
            clm.coeffs[1, l, order] = sh_coeffs[index - order]

    return clm


def grid_updater(radius, radius_arr, sh_coeffs):
    """
    Returns an initialised PySHTOOLS object given a 1D NumPy array of spherical harmonic coefficients and a radius.

    :param  int radius: The radius at which the spherical harmonic coefficients are evaluated.
    :param  np.ndarray radius_arr: 1D np array The array of radii at which the spherical harmonic coefficients are evaluated.
    :param np.ndarray sh_coeffs: 1D np array of S20RTS spherical harmonic coefficients.
    :return np.ndarray grid_data: 2D np array of shear wave velocity perturbation at the given radius.
    """

    # Find the index of the closest radius to the given radius.
    depth_idx = (np.abs(radius_arr - radius)).argmin()
    # Extract the spherical harmonic coefficients at the given radius.
    sh_coeffs_at_radius = sh_coeffs[depth_idx, :]

    # Initialise a pyshtools SHCoeffs class.
    clm = pysh.SHCoeffs.from_zeros(lmax=l_max, kind='complex', normalization='ortho')
    # Use the Condon-Shortley phase convention.
    clm.cphase = -1
    # Fill the initialised class with the spherical harmonic coefficients at the given radius.
    clm = S20RTS_to_pyshtools(sh_coeffs=sh_coeffs_at_radius, clm=clm)

    # Expand the spherical harmonic coefficients to a 2D grid.
    grid = clm.expand()

    # Return the grid values as a 2D np array.
    return grid.data


def preprocess(file_path):
    """
    Preprocesses the S20RTS spherical harmonic coefficients given a file path.

    :param str file_path: The path to the S20RTS spherical harmonic coefficients.
    :return np.ndarray sh_coeffs: 2D np array of S20RTS spherical harmonic coefficients.
    :return np.ndarray radius: 1D np array of radii at which the spherical harmonic coefficients are evaluated.
    :return np.ndarray ref_density: 2D np array of reference density values at which the spherical harmonic coefficients are evaluated.
    :return np.ndarray ref_shear_velocity: 2D np array of reference shear wave velocity values at which the spherical harmonic coefficients are evaluated.

    """

    # Load in the S20RTS spherical harmonic coefficients.
    S20RTS_data = np.loadtxt(file_path)

    radius = S20RTS_data[:, 0]
    ref_density = S20RTS_data[:, 1]
    ref_shear_velocity = S20RTS_data[:, 2]
    sh_coeffs = S20RTS_data[:, 3:]

    size = sh_coeffs.shape[1]

    # First half of the sh_coeffs array is the real part, second half is the imaginary part.
    real_sh_coeffs = sh_coeffs[:, 0: int(size / 2)]
    imag_sh_coeffs = sh_coeffs[:, int(size / 2):]

    # Convert into complex numbers.
    sh_coeffs = real_sh_coeffs + 1j * imag_sh_coeffs

    return sh_coeffs, radius, ref_density, ref_shear_velocity


l_max = 20

# Read in spherical harmonic coefficients
file_path = '../out/output.txt'
sh_coeffs, radius, ref_density, ref_shear_velocity = preprocess(file_path)

# Initial grid data evaluated at the surface.
grid_data_init = grid_updater(radius_arr=radius, radius=radius[-1], sh_coeffs=sh_coeffs)

fig, ax = plt.subplots(1, 1, figsize=(10, 5))

# Plot the data
im = ax.imshow(grid_data_init.real, cmap='viridis', extent=[-180, 180, -90, 90])
ax.set_title(fr'Re$(\delta v_s / v_s)$ at r = {radius[-1]:.2} km')

ax.set_ylabel(r'Latitude ($^\circ$)')
ax.set_xlabel(r'Longitude ($^\circ$)')

# Add colorbar for the first subplot
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r'$\delta v_s / v_s$')

# Adjust layout to prevent clipping of colorbar labels
plt.tight_layout()

# Add slider for radius
radius_ax = fig.add_axes([0.1, 0., 0.65, 0.03])
radius_slider = Slider(radius_ax, 'Radius', radius[0], radius[-1], valinit=radius[-1])


def update(val):
    # Update the shear wave velocity to be evaluated at the desired radius.
    grid_data = grid_updater(radius_arr=radius, radius=radius_slider.val, sh_coeffs=sh_coeffs)
    ax.imshow(grid_data.real, cmap='viridis', extent=[-180, 180, -90, 90])
    ax.set_title(fr'Re$(\delta v_s / v_s)$ at r = {radius_slider.val:.2} km')
    fig.canvas.draw()


radius_slider.on_changed(update)

# Show the plots
plt.show()

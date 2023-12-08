# Module to read in kernel and S20RTS data.

import numpy as np
from typing import Tuple
from validation import assert_shape

def leading_zeros(degree: int) -> str:
    """
    Adds leading zeros to the degree to match the kernel file format.

    Example:
        1  -> 001
        30 -> 030

    Parameters
    ----------
    degree : int
        The degree of the spherical harmonic.

    Returns
    -------
    str
        The degree with leading zeros.
    """
    return str(degree).zfill(3)


def parse_file(path: str, number_of_columns: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parses a file with two columns of data.

    Parameters
    ----------
    path : str
        The path to the file.
    number_of_columns : int, optional
        The number of columns in the file.

    Returns
    -------
    np.ndarray
        The first column of the file.
    np.ndarray
        The second column of the file.
    """

    with open(path, "r") as file:
        # Initialize empty lists to store the data
        radial_data = []
        kernel_data = []

        # Read and process each line
        for line in file:
            columns = line.strip().split()

            # Check if there are two columns
            if len(columns) == 2:
                try:
                    # Convert the values to floats and append to respective lists
                    radial_data.append(float(columns[0]))
                    kernel_data.append(float(columns[1]))
                except ValueError:
                    print("Error: Invalid data format on line:", line)
            else:
                print("Error: Invalid data format on line:", line)

    return np.array(radial_data), np.flip(np.array(kernel_data))

def read_in_S20RTS(file_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    This function reads in the S20RTS spherical harmonic coefficients given a file path.

    Parameters
    ----------
    file_path : str
        The path to the S20RTS spherical harmonic coefficients.

    Returns
    -------
    np.ndarray
        2D np array of S20RTS spherical harmonic coefficients. Format is

        c_{0, 0}, c_{1, -1}, c_{1, 0}, c_{1, 1}, c_{2, -2}, c_{2, -1}, ..., c_{l_max, l_max}

            .         .          .         .         .          .      ...         .
            .         .          .         .         .          .      ...         .
            .         .          .         .         .          .      ...         .
        c_{0, 0}, c_{1, -1}, c_{1, 0}, c_{1, 1}, c_{2, -2}, c_{2, -1}, ..., c_{l_max, l_max}

    np.ndarray
        1D np array of radii values.
    np.ndarray
        2D np array of reference density values in spherical harmonic format.
    np.ndarray
        2D np array of reference shear wave velocity values in spherical harmonic format.
    """

    # Load in the S20RTS spherical harmonic coefficients.
    l_max = 20

    S20RTS_data = np.loadtxt(file_path)
    if assert_shape(S20RTS_data, (None, 2 * (l_max + 1) ** 2 + 3)):
        raise ValueError(f"S20RTS_data has should have {2 * (l_max + 1) ** 2 + 3} columns, "
                         f"but it has {S20RTS_data.shape[1]} columns.")
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



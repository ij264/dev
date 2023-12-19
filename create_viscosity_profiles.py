import numpy as np
import constants
import matplotlib.pyplot as plt


def generate_viscosity_profile(depths: np.ndarray,
                               starting_viscosity: float,
                               changes_in_viscosity: np.ndarray) -> np.ndarray:
    """
    Generate a viscosity profile with a number of parameters.

    Parameters
    ----------
    depths : np.ndarray((n, ), dtype=float)
        The depths at which the viscosity changes.
    starting_viscosity : float
        The logarithm of viscosity at the surface.
    changes_in_viscosity : np.ndarray((n, ), dtype=float)
        The logarithmic change in viscosity at each depth.

    Returns
    -------
    viscosity_profile : np.ndarray
        The viscosity profile.
    """

    viscosity = np.ones(constants.number_of_radial_points) * starting_viscosity

    depth_indices = np.rint(depths / (constants.surface_boundary - constants.CMB_boundary) * constants.number_of_radial_points).astype(int)

    # Check that the number of depths is equal to the number of changes in viscosity
    assert (len(depths) == len(changes_in_viscosity))

    for i, (depth_index, change) in enumerate(zip(depth_indices, changes_in_viscosity)):
        viscosity[depth_index:] += change
    return np.power(10, viscosity)


if __name__ == "__main__":
    depths = np.array([670, ])
    delta_log_visc = np.linspace(0, 3, 2)

    for visc in delta_log_visc:
        radial_points = np.linspace(constants.surface_boundary, constants.CMB_boundary, constants.number_of_radial_points)
        viscosity = generate_viscosity_profile(depths=depths, starting_viscosity=0, changes_in_viscosity=delta_log_visc)
        plt.plot(viscosity, radial_points)
        plt.xscale('log')
        plt.show()
        file_path = f'lyness/VISC_INPUTS/step_changes/{visc:.2f}.vis'

        with open(file_path, 'w') as file:
            for i, (radius, visc) in enumerate(zip(radial_points, viscosity)):
                file.write(f'{visc}\n')
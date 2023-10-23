import numpy as np

# viscosity_depth_input = "(1, 6370), (10, 5000), (100, 3400)"

# TODO
# make sure that viscosity values are changing except for start and end.

def get_viscosities_and_depths(viscosity_depth_string: str):
    """
    Parses input string in the form '(viscosity, depth), (viscosity, depth), ...'.

    Returns two lists of floats, one for the viscosities and one for the depths.

    Parameters: viscosity_depth_string: str
                                        Input string in the form '(viscosity, depth), (viscosity, depth), ...'

    Returns: viscosities: list[float]
                          List of viscosities.

             depths: list[float]
                    List of depths.
    """

    viscosity_depth_input = viscosity_depth_string.replace('(', '').replace(')', '').split(', ')
    viscosities = [float(viscosity) for viscosity in viscosity_depth_input[::2]]
    depths = [float(depth_input) for depth_input in viscosity_depth_input[1::2]]

    return viscosities, depths


def viscosity_interpolate(viscosities, depths, root, filename):
    """
    Linear interpolation for viscosity and depth given points

    Writes a .vis file with the interpolated viscosities.

    Parameters: viscosities : array_like
                              The viscosity coordinates.
                depths : array_like
                            The depth coordinates.
                root : str
                       The root directory to write the .vis file to.
                filename : str
                           The name of the .vis file to write.

    Returns: None

    Raises: ValueError
            If the start and end depths are not 6370 km and 3400 km respectively.
    """

    if depths[-1] < depths[0]:
        depths = np.flip(depths)
        viscosities = np.flip(viscosities)

    # Ensure the start and end depths are within the range of provided data.
    if depths[0] != 3400 or depths[-1] != 6370:
        raise ValueError("Start and end depths must be 6370 km and 3400 km respectively.")

    cmb_depth = depths[0]
    surface = depths[-1]
    num_of_depth_points = 257

    # Perform linear interpolation to calculate viscosity at intermediate depths.
    interpolated_depths = np.linspace(cmb_depth, surface, num_of_depth_points)
    interpolated_viscosities = np.interp(interpolated_depths, depths, viscosities)

    # Flip the interpolated viscosities.
    interpolated_viscosities = np.flip(interpolated_viscosities)

    with open(root + filename, "w") as file:
        for viscosity in interpolated_viscosities:
            file.write(str(viscosity) + "\n")
    print(f"File written to {root + filename}.")


viscosity_depth_input = input("Enter viscosity points in the form (viscosity, depth) separated by commas: ")
viscosities, depths = get_viscosities_and_depths(viscosity_depth_input)
viscosity_interpolate(viscosities, depths, root='proto_visc/', filename='const_visc.vis')

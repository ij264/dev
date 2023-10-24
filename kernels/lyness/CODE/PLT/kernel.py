import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def get_viscosity(filepath: str):
    with open(filepath, "r") as file:
        # Initialize empty lists to store the data

        # Read and process each line
        visc_data = [float(line) for line in file]
    # convert to 2 dimensional numpy array where rows are kernel data and columns are radial data
    return np.flip(np.array(visc_data))


def save_arr(arr, file_name):
    # Saves kernels to text file
    np.savetxt(file_name, arr, delimiter=',')


def leading_zeros(degree):
    return str(degree).zfill(3)


def find_file_path(dir_to_find, root, degree, graph_type):
    file_name = graph_type + leading_zeros(degree)
    for path, dirs, files in os.walk(root):
        if dir_to_find in dirs:
            return os.path.join(path, dir_to_find + '/KERNELS/' + file_name)


def parse_file(path, number_of_columns=2):

    with open(path, "r") as file:
        # Initialize empty lists to store the data
        radial_data = []
        kernel_data = []

        # Read and process each line
        for line in file:
            # Split the line into columns using space as the delimiter
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
    # convert to 2 dimensional numpy array where rows are kernel data and columns are radial data
    return np.array(radial_data), np.flip(np.array(kernel_data))


class Kernel:
    """
    Functionality:
        - Take in a parameter which specifies whether to plot:
            - dynamic surface topography: surftopo
            - dynamic CMB topography: cmbtopo
            - geoid: geoid
            - free-air gravity anomaly: gravity
            - surface plate velocity: surfvel
            - gravity/dynamic surface topography: admittance
    """

    def __init__(self, target_dir, graph_types, degrees, root='../../', title=None):
        # TODO
        self.graph_types_dict = {
            'surftopo': 'Dynamic Surface Topography Kernel',
            'cmbtopo': 'Dynamic CMB Topography Kernel',
            'geoid': 'Geoid',
            'gravity': 'Free-air Gravity Anomaly Kernel',
            'surfvel': 'Surface Plate Velocity Kernel',
            'admittance': 'Gravity/Dynamic Surface Topography Kernel'
        }

        if not all(graph_type in self.graph_types_dict for graph_type in graph_types):
            raise ValueError(f'Invalid graph type. Please choose from the following: {self.graph_types_dict.keys()}')

        self.graph_types = graph_types
        if degrees == 'all':
            self.degrees = np.arange(1, 31)
        elif isinstance(degrees, (list, tuple, np.ndarray)):
            self.degrees = degrees

        if not all(degree in self.degrees for degree in self.degrees):
            raise ValueError(f'Only degrees 1 - 30 are supported currently.')

        self.radial_data_length = 257
        self.target_dir = target_dir
        self.title = title
        self.fig, self.axes = plt.subplots(1, len(self.graph_types) + 1, figsize=(int(10 * len(self.graph_types)), 10))
        self.kernels = np.zeros((len(self.graph_types), self.radial_data_length, len(self.degrees)))
        self.root = root
        self.kernels_to_2d_numpy_array()
        self.viscosity = get_viscosity(os.path.join(root + 'VISC_INPUTS/', 'colli.vis'))

    def kernels_to_2d_numpy_array(self):
        for i, graph_type in enumerate(self.graph_types):
            for j, degree in enumerate(self.degrees):
                path = find_file_path(self.target_dir, self.root, degree, graph_type)
                self.radial_data, self.kernel_data = parse_file(path)
                self.kernel_data = self.kernel_data
                self.kernels[i, :, j] = self.kernel_data

    def admittance(self):
        pass

    def plot(self):
        fig, axis = plt.subplots(figsize=(10, 10))
        for i, degree in enumerate(self.degrees):
            plt.plot(self.kernels[0, :, i], self.radial_data, label=f'l = {degree}')
        axis.xaxis.tick_top()
        axis.xaxis.set_label_position('top')
        axis.set_xlabel(self.graph_types[0])
        axis.set_ylabel('Radius (km)')
        axis.legend()
        plt.show()

    def imshow(self):
        self.axes[0].plot(self.viscosity, self.radial_data, color='black')
        self.axes[0].set_xscale('log')
        self.axes[0].set_xlabel('$\mu \,\, (Pa \, s)$')
        self.axes[0].set_ylabel('Radius (km)')
        self.axes[0].set_ylim(self.radial_data.min(), self.radial_data.max())
        cmap = LinearSegmentedColormap.from_list('CustomColors', ['blue', 'white', 'red'], N=256)
        for i, graph_type in enumerate(self.graph_types):
            # cmap = LinearSegmentedColormap.from_list('CustomColors', ['blue', 'white'], N=256)
            im = self.axes[i + 1].imshow(self.kernels[i, :, :],
                                         extent=[np.min(self.degrees), np.max(self.degrees), np.min(self.radial_data),
                                                 np.max(self.radial_data)], aspect='auto', cmap=cmap)
            self.axes[i].set_ylim(self.radial_data.min(), self.radial_data.max())
            cbar = self.fig.colorbar(im, ax=self.axes[i + 1], shrink=0.7)
            self.axes[i + 1].set_xlabel('Degree')
            self.axes[i + 1].set_ylabel('Radius (km)')
            self.axes[i + 1].set_title(self.graph_types[i])
        plt.tight_layout()
        plt.show()

    def contour(self):

        # Create a regular grid to interpolate the data.
        x = np.arange(1, 31)
        y = np.arange(1, 258)
        X, Y = np.meshgrid(x, y)

        im = self.ax.imshow(self.kernels, cmap='red', extent=[0, 30, 0, 258], aspect='auto')
        CS = self.ax.contour(X, Y, self.kernels)
        self.ax.clabel(CS, inline=True, fontsize=10)
        # Add a colorbar for reference
        # plt.colorbar()

        plt.xlabel('X Axis')
        plt.ylabel('Y Axis')
        plt.title('Interpolated Contour Plot')

        plt.show()
        plt.clabel(CS, fontsize=10, colors='k')
        self.ax.set_xlabel('Degree')
        self.ax.set_ylabel('Radius (km)')
        self.ax.set_title(self.title)
        plt.show()

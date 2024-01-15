import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize

def get_viscosity(filepath: str):
    with open(filepath, "r") as file:
        # Initialize empty lists to store the tomographic_models

        # Read and process each line
        visc_data = [float(line) for line in file]
    # convert to 2 dimensional numpy array where rows are kernel tomographic_models and columns are radial tomographic_models
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
        # Initialize empty lists to store the tomographic_models
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
                    print("Error: Invalid tomographic_models format on line:", line)
            else:
                print("Error: Invalid tomographic_models format on line:", line)
    # convert to 2 dimensional numpy array where rows are kernel tomographic_models and columns are radial tomographic_models
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
        # self.fig, self.axes = plt.subplots(1, len(self.graph_types) + 1, figsize=(int(8 * len(self.graph_types)), 10))
        self.cmap = LinearSegmentedColormap.from_list('CustomColors', ['blue', 'white', 'red'], N=256)
        self.kernels = np.zeros((len(self.graph_types), self.radial_data_length, len(self.degrees)))
        self.root = root
        self.radial_data, self.kernels = self.kernels_to_2d_numpy_array(self.kernels, self.graph_types, self.degrees)
        self.viscosity = get_viscosity(os.path.join(root + self.target_dir, 'viscosity_profile.vis'))

    def kernels_to_2d_numpy_array(self, kernels, graph_types, degrees):
        for i, graph_type in enumerate(graph_types):
            for j, degree in enumerate(degrees):
                path = find_file_path(self.target_dir, self.root, degree, graph_type)
                radial_data, kernel_data = parse_file(path)
                kernels[i, :, j] = kernel_data
        return radial_data, kernels

    def manual_admittance(self):
        kernels = np.zeros((2, self.radial_data_length, len(self.degrees)))
        _, gravity_and_surftopo = self.kernels_to_2d_numpy_array(kernels, graph_types=['gravity', 'surftopo'],
                                                                 degrees=self.degrees)
        gravity_kernel = gravity_and_surftopo[0, :, :]
        surftopo_kernel = gravity_and_surftopo[1, :, :]
        admittance_kernel = gravity_kernel / surftopo_kernel
        fig, ax = plt.subplots(figsize=(8, 8))
        return admittance_kernel

    def plot(self):
        match self.graph_types[0]:
            case 'surftopo':
                lims = [-1.1, 0.1]
            case 'geoid':
                lims = [-0.5, 0.5]
        linestyles = ['-', '--', ':']
        fig, axes = plt.subplots(1, 2, figsize=(6, 4))
        axes[0].plot(self.viscosity, self.radial_data, color='black')
        axes[0].set_xscale('log')
        axes[0].set_xlabel('$\eta/\eta_0$')
        axes[0].set_ylabel('Radius (km)')
        axes[0].set_ylim(self.radial_data.min(), self.radial_data.max())
        axes[0].set_xlim([0.5, 5e4])

        for i, (linestyle, degree) in enumerate(zip(linestyles, self.degrees)):
            axes[1].plot(np.flip(self.kernels[0, :, i]), self.radial_data, label=f'l = {degree}', linestyle=linestyle, color='black')
            axes[1].set_xlim(lims)

        for i, ax in enumerate(fig.axes):
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position('top')
            ax.set_ylim(self.radial_data.min(), self.radial_data.max())

            ax.set_ylabel('Radius (km)')
            ax.set_xlabel(self.graph_types[0] + ' Kernel')
            #
            # axes[i].invert_xaxis()

        plt.legend()

        plt.tight_layout()
        plt.savefig(f'figures/{self.target_dir[7:]}_{self.graph_types[0]}.png', dpi=300)
        plt.show()

    def imshow(self, filename, normalise=True, show_viscosity=True, save=False):
        if self.kernels.ndim == 2:
            kernels = np.expand_dims(self.kernels, axis=0)

        number_of_subplots = self.kernels.shape[0]
        if show_viscosity:
            number_of_subplots += 1
        columns = int(np.ceil(np.sqrt(number_of_subplots)))
        rows = int(np.ceil(number_of_subplots / columns))
        fig, axes = plt.subplots(rows, columns, figsize=(10, 10), squeeze=0)

        for i, ax in enumerate(fig.axes):
            if i > number_of_subplots - 1:
                ax.set_visible(False)

        for i, ax in enumerate(fig.axes):
            if i == 0 and show_viscosity:
                ax.plot(self.viscosity, self.radial_data, color='black')
                ax.set_xscale('log')
                ax.set_xlabel('$\mu \,\, (Pa \, s)$')
                ax.set_ylabel('Radius (km)')
                ax.set_ylim(self.radial_data.min(), self.radial_data.max())
            elif i > number_of_subplots - 1:
                break
            else:
                if normalise:
                    vmin = -1
                    vmax = 1

                else:
                    vmin = np.min(self.kernels[i - 1, :, :])
                    vmax = np.max(self.kernels[i - 1, :, :])

                norm = Normalize(vmin, vmax)

                im = ax.imshow(self.kernels[i-1, :, :],
                               extent=[np.min(self.degrees), np.max(self.degrees), np.min(self.radial_data),
                                       np.max(self.radial_data)], aspect='auto', cmap='seismic', norm=norm)
                ax.set_ylim(self.radial_data.min(), self.radial_data.max())
                cbar = self.fig.colorbar(im, ax=ax, shrink=0.7)
                ax.set_xlabel('Degree')
                ax.set_title(self.graph_types[i-1])
                ax.set_ylabel('Radius (km)')
        plt.tight_layout()
        if save:
            plt.savefig(f'{filename}.png', dpi=300)
        plt.show()

    def surface_admittance(self):
        admittance = np.zeros((1, 257, 30))
        _, admittance = self.kernels_to_2d_numpy_array(admittance, ['admittance'], [i for i in range(1, 31)])
        return admittance[0, 0, :]


    def contour(self):

        # Create a regular grid to interpolate the tomographic_models.
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

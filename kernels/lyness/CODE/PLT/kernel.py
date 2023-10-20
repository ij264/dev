import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.interpolate
from matplotlib.colors import LinearSegmentedColormap


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

        # try:
        #     if self.graph_type == 'surftopo':
        #         self.xlabel = 'Dynamic Surface Topography Kernel'
        #     elif self.graph_type == 'cmbtopo':
        #         self.xlabel = 'Dynamic CMB Topography Kernel'
        #     elif self.graph_type == 'geoid':
        #         self.xlabel = 'Geoid'
        #     elif self.graph_type == 'gravity':
        #         self.xlabel = 'Free-air Gravity Anomaly Kernel'
        #     elif self.graph_type == 'surfvel':
        #         self.xlabel = 'Surface Plate Velocity Kernel'
        #     elif self.graph_type == 'admittance':
        #         self.xlabel = 'Gravity/Dynamic Surface Topography Kernel'
        # except NameError:
        #     print(graph_type, 'is not a valid kernel here.')
        self.visc_profile = self.parse_file(root + 'INPUT/visc_profile.txt')[1]
        self.degrees = degrees
        self.radial_data_length = 257
        self.target_dir = target_dir
        self.title = title
        self.fig, self.axes = plt.subplots(1, len(self.graph_types), figsize=(10, 10))
        self.kernels = np.zeros((len(self.graph_types), self.radial_data_length, len(self.degrees)))
        self.root = root
        self.kernels_to_2d_numpy_array()

    def kernels_to_2d_numpy_array(self):
        for i, graph_type in enumerate(self.graph_types):
            for j, degree in enumerate(self.degrees):
                path = self.find_file_path(self.target_dir, self.root, degree, graph_type)
                self.radial_data, self.kernel_data = self.parse_file(path)
                self.kernels[i, :, j] = self.kernel_data
        # Return the kernel data as a 2D numpy array of size (degrees, radial_data)
        # self.kernels = np.transpose(np.array(self.kernels))
        # self.kernels = np.flip(self.kernels, axis=0)

    def plot(self, style='line', marker_size=5):
        # Initialise figure and axes
        # print self.degrees
        for i, degree in enumerate(self.degrees):
            if style == 'line':
                self.ax.plot(self.kernels[:, i], self.radial_data, label=f'l = {degree}')
            elif style == 'scatter':
                self.ax.scatter(self.kernels[:, i], self.radial_data, label=f'l = {degree}', s=marker_size)
        self.ax.xaxis.tick_top()
        self.ax.xaxis.set_label_position('top')
        self.graph_type = self.graph_type
        self.ylabel = 'Radius (km)'

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)
        self.ax.legend()
        plt.show()

    def imshow(self):
        cmap = LinearSegmentedColormap.from_list('CustomColors', ['blue', 'white', 'red'], N=256)
        for i, graph_type in enumerate(self.graph_types):

            # cmap = LinearSegmentedColormap.from_list('CustomColors', ['blue', 'white'], N=256)
            im = self.axes[i].imshow(self.kernels[i, :, :], extent=[np.min(self.degrees), np.max(self.degrees), np.min(self.radial_data),
                                                      np.max(self.radial_data), ], aspect='auto', cmap=cmap)
            cbar = self.fig.colorbar(im, ax=self.axes[i], shrink=0.7)
            self.axes[i].set_xlabel('Degree')
            self.axes[i].set_ylabel('Radius (km)')
            self.axes[i].set_title(self.title)
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

    def save(self, arr, file_name):
        # Saves kernels to text file
        np.savetxt(file_name, arr, delimiter=',')

    def show(self):
        pass

    def leading_zeros(self, degree):
        return str(degree).zfill(3)

    def find_file_path(self, dir_to_find, root, degree, graph_type):
        file_name = graph_type + self.leading_zeros(degree)
        for path, dirs, files in os.walk(root):
            if dir_to_find in dirs:
                return os.path.join(path, dir_to_find + '/KERNELS/' + file_name)

    def parse_file(self, path):

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
        return np.array(radial_data), np.array(kernel_data)

import matplotlib.pyplot as plt
import os


class Graph:
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

    def __init__(self, dir_to_find, graph_type, degrees, root='../../', title=None):
        # TODO

        self.graph_type = graph_type
        self.ylabel = 'Radius (km)'
        self.xlabel = None

        if self.graph_type == 'surftopo':
            self.xlabel = 'Dynamic Surface Topography Kernel'
        elif self.graph_type == 'cmbtopo':
            self.xlabel = 'Dynamic CMB Topography Kernel'
        elif self.graph_type == 'geoid':
            self.xlabel = 'Geoid'
        elif self.graph_type == 'gravity':
            self.xlabel = 'Free-air Gravity Anomaly Kernel'
        elif self.graph_type == 'surfvel':
            self.xlabel = 'Surface Plate Velocity Kernel'
        elif self.graph_type == 'admittance':
            self.xlabel = 'Gravity/Dynamic Surface Topography Kernel'
        else:
            print(graph_type, 'is not a valid kernel here.')

        # Initialise figure and axes
        self.fig, self.ax = plt.subplots(figsize=(8, 10))
        self.ax.xaxis.tick_top()
        self.ax.xaxis.set_label_position('top')

        for degree in degrees:
            self.degree = degree
            self.path = self.find_file_path(dir_to_find, root, self.degree, self.graph_type)
            self.radial_data, self.kernel_data = self.parse_files(self.path)
            self.plot(self.kernel_data, self.radial_data)

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(title)
        self.ax.legend()
        plt.show()

    def plot(self, x, y, style='line', marker_size=5):
        if style == 'line':
            self.ax.plot(x, y, label=f'l = {self.degree}')
        elif style == 'scatter':
            self.ax.scatter(x, y, label=f'l = {self.degree}', s=marker_size)

    def show(self):
        pass

    def leading_zeros(self, degree):
        return str(degree).zfill(3)

    def find_file_path(self, dir_to_find, root, degree, graph_type):
        file_name = graph_type + self.leading_zeros(degree)
        for path, dirs, files in os.walk(root):
            if dir_to_find in dirs:
                return os.path.join(path, dir_to_find + '/KERNELS/' + file_name)

    def parse_files(self, path):

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

        return radial_data, kernel_data

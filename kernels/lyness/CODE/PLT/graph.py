import matplotlib.pyplot as plt


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

    def __init__(self, graph_type, x_data, y_data, degree, title=None):

        self.graph_type = graph_type
        self.ylabel = 'Radius (km)'
        if self.graph_type == 'surftopo':
            self.xlabel = 'Dynamic Surface Topography'
        elif self.graph_type == 'cmbtopo':
            self.xlabel = 'Dynamic CMB Topography'
        elif self.graph_type == 'geoid':
            self.xlabel = 'Geoid'
        elif self.graph_type == 'gravity':
            self.xlabel = 'Free-air Gravity Anomaly'
        elif self.graph_type == 'surfvel':
            self.xlabel = 'Surface Plate Velocity'
        elif self.graph_type == 'admittance':
            self.xlabel = 'Gravity/Dynamic Surface Topography'

        # Initialise figure and axes
        self.fig, self.ax = plt.subplots()
        self.x = x_data
        self.y = y_data
        self.degree = degree
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(title)
        self.ax.legend()

    def plot(self, style='line', label=None, marker_size=5):
        if type == 'line':
            self.ax.plot(self.x, self.y, label=label)
        elif type == 'scatter':
            self.ax.scatter(self.x, self.y, label=label, s=marker_size)

    def show(self):
        plt.show()

import os
import numpy as np
import xarray as xr

def leading_zeros(degree):
    return str(degree).zfill(3)


def find_kernel_path(dir_to_find, root, degree, graph_type):
    file_name = graph_type + leading_zeros(degree)
    for path, dirs, files in os.walk(root):
        if dir_to_find in dirs:
            return os.path.join(path, dir_to_find + '/KERNELS/' + file_name)


def parse_file(path, filetype, number_of_columns=2):
    if filetype == 'np':
        return np.loadtxt(path)
    elif filetype == 'nc':
        ds = xr.open_dataarray(path)
        return ds.to_dataframe()
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


def kernel_at_degree(path, degree):
    radial_data, kernel_data = parse_file(path)
    return kernel_data

from kernel import Kernel
root = 'lyness/'
target_dir = 'OUTPUT_const_visc'
graph_types = ['surftopo']
degrees = [2, 5, 10]


Ker = Kernel(target_dir=target_dir, graph_types=graph_types, degrees=degrees, root=root)
Ker.plot()
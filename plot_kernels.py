from kernel import Kernel
root = 'lyness/'
target_dir = 'OUTPUT_F10a'
graph_types = ['surftopo']
degrees = [2, 20, 30, 40]


Ker = Kernel(target_dir=target_dir, graph_types=graph_types, degrees=degrees, root=root)
Ker.plot()
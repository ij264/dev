import numpy as np

a = np.array([[[1, 1], [2, 2]], [[3, 3], [4, 4]]])
b = np.array([0.5, 2])

print(a * b[:, None, None])
# for i, _ in enumerate(b):

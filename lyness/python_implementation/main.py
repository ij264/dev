import numpy as np
from  scipy import linalg


def construct_A_matrix(l: int):

    L = l * (l + 1)

    A = np.array(
        [
            [-2, L, 0, 0, 0, 0],
            [-1, 1, 0, 1, 0, 0],
            [12, -6*L, 1, L, 0, -1],
            [-6, 2*(2*L-1), -1, -2, -1, 0],
            [0, 0, 0, 0, 1, 1],
            [0, 0, 0, -0, L, 1]
        ]
    )
    return A

def exponentiate(A):
    return linalg.expm(A)

A = construct_A_matrix(2)

print(exponentiate(A))
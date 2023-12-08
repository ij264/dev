from unittest import TestCase
import unittest
import numpy as np
import process
import random

def generate_test_coeffs(num_radial, l_max):
    """
    Returns a 3D np array of test coefficients of the form

    [
    [00 00 00 ... ... ... ... 00]
    [10 11 00 ... ... ... ... 00]
    [20 21 22 ... ... ... ... 00]
    [30 31 32 ... ... ... ... 00]
    ...
    [l0 ... ... ... .. ... ... ll],

    [00 00 00 ... ... ... ... 00],
    [00 -10 -11 ... ... ... ... 00],
    [00 -20 -21 ... ... ... ... 00],
    ...
    [00 -l0 ... ... ... .. ... ... -ll]
    ] x number of radial points.
    """

    # radial points x 2 x l_max+1 x l_max+1
    mat = np.array([[[[0 for i in range(l_max+1)] for j in range(l_max+1)] for k in range(2)]
                     for r in (range(num_radial))])

    # Positive orders.
    for r in range(num_radial):
        for l in range(l_max+1):
            for m in range(l+1):
                mat[r][0][l][m] = int(str(l) + str(m))

    # Negative orders.
    for r in range(num_radial):
        for l in range(l_max+1):
            for m in range(1, l+1):
                mat[r][1][l][m] = -int(str(l) + str(m))


    return mat

class TestIntegrate(TestCase):
    def setUp(self):
        self.radius_min = 3.480e3
        self.radius_max = 6.470e3
        self.radial_points = 185
        self.radius = np.linspace(self.radius_min,
                                  self.radius_max,
                                  self.radial_points)
        self.l_max = 20
        self.test_coeffs = generate_test_coeffs(self.radial_points,
                                                     l_max=self.l_max)
        self.test_kernel = np.ones_like(self.radius)


class TestInit(TestIntegrate):
    def test_integration(self):
        degree = random.randint(1, self.l_max)
        integral = process.integrate(l=degree,
                       density_sh_lm=self.test_coeffs,
                       kernel=self.test_kernel,
                       radius_arr=self.radius)
        print('Done!')

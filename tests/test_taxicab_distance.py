# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 23:19:28 2022

@author: Cal
"""

import os
os.chdir("..")

import unittest
import math

import leod


class TestTaxicab(unittest.TestCase):

    def test_taxicab_sphere(self):
        print("Test spherical taxicab distance")
        conv = math.pi / 180.0
        start_point = [50.0*conv, 60.0*conv]
        end_point = [90.0*conv, 0.0*conv]
        d = leod.geo.taxicab.sphere_td(1.0, start_point, end_point)
        self.assertAlmostEqual(d, 1.500331566, 9, "Expect d = 1.500331566")
        print("Test passed")
        print(" ")







if __name__ == '__main__':
    unittest.main()
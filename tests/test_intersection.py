# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 14:09:52 2022

@author: Cal
"""

import os
os.chdir("..")

import unittest
import math
import numpy as np

import leod

class TestIntersection(unittest.TestCase):

    def test_ellipsoid_plane(self):
        print("Test ellipsoid-plane intersection")
        shape = leod.shape.EllipsoidShape(1.0, 1.0, 1.0)
        n = np.array([0.0, 0.0, 1.0])
        q = np.array([0.0, 0.0, 0.0])
        c, a, b = leod.intersection.ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 1.0, 9, "Expect a_axis = 1")
        self.assertAlmostEqual(b, 1.0, 9, "Expect b_axis = 1")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        n = np.array([1.0, 0.0, 0.0])
        q = np.array([0.0, 0.0, 0.0])
        c, a, b = leod.intersection.ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 1.0, 9, "Expect a_axis = 1")
        self.assertAlmostEqual(b, 1.0, 9, "Expect b_axis = 1")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        n = np.array([0.0, 2.0, 0.0])
        q = np.array([0.0, 0.0, 0.0])
        c, a, b = leod.intersection.ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 1.0, 9, "Expect a_axis = 1")
        self.assertAlmostEqual(b, 1.0, 9, "Expect b_axis = 1")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        shape = leod.shape.EllipsoidShape(2.0, 2.0, 1.0)
        n = np.array([1.0, 0.0, 0.0])
        q = np.array([0.0, 0.0, 0.0])
        c, a, b = leod.intersection.ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 1.0, 9, "Expect a_axis = 1")
        self.assertAlmostEqual(b, 2.0, 9, "Expect b_axis = 2")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        n = np.array([0.0, 0.0, 1.0])
        q = np.array([0.0, 0.0, 0.5])
        c, a, b = leod.intersection.ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, math.sqrt(3.0), 9, "Expect a_axis = sqrt(3)")
        self.assertAlmostEqual(b, math.sqrt(3.0), 9, "Expect b_axis = sqrt(3)")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        print("Test passed")
        print(" ")







if __name__ == '__main__':
    unittest.main()
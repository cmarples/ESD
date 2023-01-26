# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 14:09:52 2022

@author: Cal
"""

import os
os.chdir("..")

import unittest
from math import sqrt
from numpy import array
from esd.shape import EllipsoidShape
from esd.intersection import ellipsoid_plane

class TestIntersection(unittest.TestCase):

    def test_ellipsoid_plane(self):
        print("Test ellipsoid-plane intersection")
        shape = EllipsoidShape(1.0, 1.0, 1.0)
        n = array([0.0, 0.0, 1.0])
        q = array([0.0, 0.0, 0.0])
        c, a, b = ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 1.0, 9, "Expect a_axis = 1")
        self.assertAlmostEqual(b, 1.0, 9, "Expect b_axis = 1")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        n = array([1.0, 0.0, 0.0])
        q = array([0.0, 0.0, 0.0])
        c, a, b = ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 1.0, 9, "Expect a_axis = 1")
        self.assertAlmostEqual(b, 1.0, 9, "Expect b_axis = 1")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        n = array([0.0, 2.0, 0.0])
        q = array([0.0, 0.0, 0.0])
        c, a, b = ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 1.0, 9, "Expect a_axis = 1")
        self.assertAlmostEqual(b, 1.0, 9, "Expect b_axis = 1")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        shape = EllipsoidShape(2.0, 2.0, 1.0)
        n = array([1.0, 0.0, 0.0])
        q = array([0.0, 0.0, 0.0])
        c, a, b = ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, 2.0, 9, "Expect a_axis = 2")
        self.assertAlmostEqual(b, 1.0, 9, "Expect b_axis = 1")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        n = array([0.0, 0.0, 1.0])
        q = array([0.0, 0.0, 0.5])
        c, a, b = ellipsoid_plane(shape, n, q)
        self.assertAlmostEqual(a, sqrt(3.0), 9, "Expect a_axis = sqrt(3)")
        self.assertAlmostEqual(b, sqrt(3.0), 9, "Expect b_axis = sqrt(3)")
        self.assertAlmostEqual(c[0], 0.0, 9, "Expect center at the origin")
        self.assertAlmostEqual(c[1], 0.0, 9, "Expect center at the origin")
        
        print("Test passed")
        print(" ")







if __name__ == '__main__':
    unittest.main()
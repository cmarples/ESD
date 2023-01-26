# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 12:54:06 2022

@author: Callum Marples

This file is used for implementing unit tests of coordinate conversions
"""

import os
os.chdir("..")

import unittest
import math

from esd.shape import EllipsoidShape

class TestConversions(unittest.TestCase):

    def test_polars(self):
        print("Test conversion between Cartesian and scaled spherical polar (theta, phi) coordinates")
        a = 3.0
        b = 2.0
        c = 1.0
        E = EllipsoidShape(a, b, c)
        # Start with polar coordinates
        th = 75.0 * math.pi / 180.0
        ph = -35.0 * math.pi / 180.0
        # Convert to Cartesians
        [x, y, z] = E.polar2cart(th, ph)
        self.assertAlmostEqual(x, 2.373720346, 9, "Expect x = 2.373720346")
        self.assertAlmostEqual(y, -1.108064586, 9, "Expect y = -1.108064586")
        self.assertAlmostEqual(z, 0.258819045, 9, "Expect x = 0.258819045")
        # Convert back to polars
        [th_out, ph_out] = E.cart2polar(x, y, z)
        self.assertAlmostEqual(th_out*180.0/math.pi, 75.0, 9, "Expect th = 75 degrees")
        # phi should be given in the range [0, 360]
        self.assertAlmostEqual(ph_out*180.0/math.pi, 325.0, 9, "Expect ph = 325 degrees")
        print("Test passed")
        print(" ")
        
    def test_ellipsoidal(self):
        print("Test conversion between Cartesian and ellipsoidal (beta, lambda) coordinates")
        a = 3.0
        b = 2.0
        c = 1.0
        E = EllipsoidShape(a, b, c)
        # Start with ellipsoidal coordinates
        be = 75.0 * math.pi / 180.0
        lm = -35.0 * math.pi / 180.0
        # Convert to Cartesians
        [x, y, z] = E.ellip2cart(be, lm)
        self.assertAlmostEqual(x, 1.981447713, 9, "Expect x = 1.981447713")
        self.assertAlmostEqual(y, -0.296905011, 9, "Expect y = -0.296905001")
        self.assertAlmostEqual(z, 0.7360194474, 9, "Expect z = 0.7360194474")
        # Convert back to polars
        [be_out, lm_out, d] = E.cart2ellip(x, y, z)
        self.assertAlmostEqual(be_out*180.0/math.pi, 75.0, 9, "Expect th = 75 degrees")
        self.assertAlmostEqual(lm_out*180.0/math.pi, -35.0, 9, "Expect th = -35 degrees")
        print("Test passed")
        print(" ")
            
if __name__ == '__main__':
    unittest.main()


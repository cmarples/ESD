# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 17:29:53 2022

@author: Callum Marples

This file is used for implementing unit tests.
"""

import unittest

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid

class TestGeo(unittest.TestCase):

    def test_geo(self):
        a = 3.0
        b = 2.0
        c = 1.0
        E = EllipsoidShape(a, b, c)
        G = GeoGrid(E, 19, 36)
        
        self.assertEqual(G.find_pixel_index(12, 19), 416, "Should be 416")
        self.assertEqual(G.find_theta_index(416), 12, "Should be 12")
        self.assertEqual(G.find_phi_index(416, 12), 19, "Should be 19")
        
        
if __name__ == '__main__':
    unittest.main()
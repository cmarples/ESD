# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 17:29:53 2022

@author: Callum Marples

This file is used for implementing unit tests.
"""

import os
os.chdir("..")

import unittest
import math

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid
from leod.geo_fmm import GeoFMM

class TestGeo(unittest.TestCase):

    def test_geo(self):
        a = 3.0
        b = 2.0
        c = 1.0
        E = EllipsoidShape(a, b, c)
        G = GeoGrid(E, 19, 36)
        
        # Test pixel finders
        self.assertEqual(G.get_pixel_index(12, 19), 416, "Should be 416")
        self.assertEqual(G.get_theta_index(416), 12, "Should be 12")
        self.assertEqual(G.get_phi_index(416, 12), 19, "Should be 19")
        
        # Test neighbour distance (using poles and symmetry)
        neighbour = 0
        d_north = G.get_distance(G.pixel[0], neighbour)
        neighbour = 35 # Symmetry with north pole
        d_south = G.get_distance(G.pixel[613], neighbour)
        self.assertAlmostEqual(d_north, d_south, 7, "By symmetry, these should be equal")
        
        # Test geodesic distances for a particular example
        grid = GeoGrid(E, 181, 360)
        th_0 = 90.0 * math.pi / 180.0
        ph_0 = 0.0
        fmm0 = GeoFMM(grid, th_0, ph_0)
        fmm1 = GeoFMM(grid, th_0, ph_0)
        fmm2 = GeoFMM(grid, th_0, ph_0)
        
        th_1 = 50.0 * math.pi / 180.0
        ph_1 = 60.0 * math.pi / 180.0
        
        d0 = fmm0.calculate_geodesics(0, th_1, ph_1)
        d1 = fmm1.calculate_geodesics(1, th_1, ph_1)
        d2 = fmm2.calculate_geodesics(2, th_1, ph_1)
        
        self.assertAlmostEqual(d0, 2.8676959922, 7, "Expect d = 2.8676959922 for 4-neighbour Dijkstra")
        self.assertAlmostEqual(d1, 2.3854095543, 7, "Expect d = 2.3854095543 for 1st order FMM")
        self.assertAlmostEqual(d2, 2.3550639686, 7, "Expect d = 2.3550639686 for 2nd order FMM")
        
if __name__ == '__main__':
    unittest.main()
    
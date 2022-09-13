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
from leod.sphere_geodesics import great_circle_distance
from leod.triaxial_geodesics import boundary_value_method
try:
    from geographiclib.geodesic import Geodesic
    from leod.spheroid_geodesics import spheroid_geo_distance
    geo_flag = True
    print(" ")
    print("GeographicLib is installed")
except ModuleNotFoundError:
    geo_flag = False
    print(" ")
    print("GeographicLib is not installed. Tests involving this library will not be run.")
print(" ")    

class TestGeo(unittest.TestCase):
    
    def test_8_neighbour(self):
        print("Simple tests of 8-neighbour GeoGrid")
        a = 3.0
        a = 3.0
        b = 2.0
        c = 1.0
        E = EllipsoidShape(a, b, c)
        G = GeoGrid(E, 5, 8, neighbour8=True)
        
        self.assertEqual(G.no_pixels, 26, "Should be 26")
        self.assertEqual(G.pixel[9].neighbour, [1, 17, 16, 10, 8, 2, 24, 18], "Should be [1, 17, 16, 10, 8, 2, 24, 18]")
        self.assertEqual(G.pixel[8].neighbour, [0, 16, 7, 1, 0, 0, 15, 9], "Should be [0, 16, 7, 1, 0, 0, 15, 9]")
        self.assertEqual(G.pixel[20].neighbour, [12, 25, 19, 21, 11, 13, 25, 25], "Should be [12, 25, 19, 21, 11, 13, 25, 25]")
        print("Test passed")
        print(" ")
            
    def test_geo(self):
        print("Simple tests for initialisation of EllipsoidShape and GeoGrid objects")
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
        
    def test_geographiclib_sphere(self):
        if geo_flag == True:
            print(" ")
            print("Test use of GeographicLib on a spherical example")
            E = EllipsoidShape(1.0, 1.0, 1.0)
            th_0 = 90.0 * math.pi / 180.0
            ph_0 = 0.0
            th_1 = 50.0 * math.pi / 180.0
            ph_1 = 60.0 * math.pi / 180.0
            
            c = great_circle_distance(1.0, th_0, ph_0, th_1, ph_1)
            s = spheroid_geo_distance(E, th_0, ph_0, th_1, ph_1)
            self.assertAlmostEqual(s, c, 7, "Expect d = 1.17773051445231")
            print("Test passed")
            print(" ")
    
    def test_geographiclib_spheroid(self):
        if geo_flag == True:
            print(" ")
            print("Test use of GeographicLib on a spheroid, comparing to a WGS84 example from the GeographicLib documentation.")
            E = EllipsoidShape(6378137.0, 6378137.0, 6356752.3142) # WGS84 ellipsoid
            th_0 = (90.0 - 41.32) * math.pi / 180.0
            ph_0 = 174.81 * math.pi / 180.0
            th_1 = (90.0 + 40.96) * math.pi / 180.0
            ph_1 = -5.50 * math.pi / 180.0
            s = spheroid_geo_distance(E, th_0, ph_0, th_1, ph_1)      
            self.assertAlmostEqual(s, 19959679.26735382, 3, "Expect d = 19959679.267 m")
            print("Test passed")
            print(" ")
            
    
            
if __name__ == '__main__':
    unittest.main()
    
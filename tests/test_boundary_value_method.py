# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:58:10 2022

@author: Callum Marples

This file is used for implementing unit tests of the boundary value method
of Panou.

The tests reproduce the values in Tables 1-4 of Panou 2013, to ensure that the
implementation used in LEOD is correct.
"""

import os
os.chdir("..")

import unittest
import math

from esd.shape import EllipsoidShape
from esd.triaxial_geodesics import boundary_value_method

print(" ")    

class TestBVM(unittest.TestCase):
    
        
    
    def test_boundary_value_method_table_1(self):
        print(" ")
        print("Test use of the boundary value method on a triaxial ellipsoid, comparing to Table 1 of Panou, 2013")
        a = 6378172.0
        b = 6378103.0 # Earth parameters
        c = 6356753.0
        E = EllipsoidShape(a, b, c)
        conv = math.pi / 180.0 # Conversion factor from degrees to radians
        dp = 4 # Number of decimal places used in comparison
        
        s = boundary_value_method(E, 0.0, 0.0, 0.0, 90.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 10018754.9569, dp, "Expect d = 10018754.9569 m")
        
        s = boundary_value_method(E, 1.0*conv, 0.0, -80.0*conv, 5.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 8947130.7221, dp, "Expect d = 8947130.7221 m")
        
        s = boundary_value_method(E, 5.0*conv, 0.0, -60.0*conv, 40.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 8004762.4330, dp, "Expect d = 8004762.4330 m")
        
        s = boundary_value_method(E, 30.0*conv, 0.0, -30.0*conv, 175.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 19547128.7971, dp, "Expect d = 19547128.7971 m")
        
        s = boundary_value_method(E, 60.0*conv, 0.0, 60.0*conv, 175.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 6705715.1610, dp, "Expect d = 6705715.1610 m")
        
        s = boundary_value_method(E, 75.0*conv, 0.0, 80.0*conv, 120.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 2482501.2608, dp, "Expect d = 2482501.2608 m")
        
        s = boundary_value_method(E, 80.0*conv, 0.0, 60.0*conv, 90.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 3519745.1283, dp, "Expect d = 3519745.1283 m")
        
        print("Test passed")
        print(" ")
    
    
    
    def test_boundary_value_method_table_2(self):
        print(" ")
        print("Test use of the boundary value method on a triaxial ellipsoid, comparing to Table 2 of Panou, 2013")
        a = 6378172.0
        b = 6378103.0 # Earth parameters
        c = 6356753.0
        E = EllipsoidShape(a, b, c)
        conv = math.pi / 180.0 # Conversion factor from degrees to radians
        dp = 4 # Number of decimal places used in comparison
        
        s = boundary_value_method(E, 0.0, -90.0*conv, 0.0, 89.5*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 19981849.8629, dp, "Expect d = 19981849.8629 m")
        
        s = boundary_value_method(E, 1.0*conv, -90.0*conv, 1.0*conv, 89.5*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 19776667.0342, dp, "Expect d = 19776667.0342 m")
        
        s = boundary_value_method(E, 5.0*conv, -90.0*conv, 5.0*conv, 89.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 18889165.0873, dp, "Expect d = 18889165.0873 m")
        
        s = boundary_value_method(E, 30.0*conv, -90.0*conv, 30.0*conv, 86.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 13331814.6078, dp, "Expect d = 13331814.6078 m")
        
        s = boundary_value_method(E, 60.0*conv, -90.0*conv, 60.0*conv, 78.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 6637321.6350, 3, "Expect d = 6637321.6350 m")
        
        s = boundary_value_method(E, 75.0*conv, -90.0*conv, 75.0*conv, 66.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 3267941.2812, dp, "Expect d = 3267941.2812 m")
        
        s = boundary_value_method(E, 80.0*conv, -90.0*conv, 80.0*conv, 55.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 2132316.9048, dp, "Expect d = 2132316.9048 m")
        
        print("Test passed")
        print(" ")
    
    
    def test_boundary_value_method_table_3(self):
        print(" ")
        print("Test use of the boundary value method on a triaxial ellipsoid, comparing to Table 3 of Panou, 2013")
        a = 6378172.0
        b = 6378103.0 # Earth parameters
        c = 6356753.0
        E = EllipsoidShape(a, b, c)
        conv = math.pi / 180.0 # Conversion factor from degrees to radians
        dp = 4 # Number of decimal places used in comparison
        
        s = boundary_value_method(E, 0.0, 0.5*conv, 80.0*conv, 0.5*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 8831874.3717, dp, "Expect d = 8831874.3717 m")
        
        s = boundary_value_method(E, -1.0*conv, 5.0*conv, 75.0*conv, 5.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 8405370.4947, dp, "Expect d = 8405370.4947 m")
        
        s = boundary_value_method(E, -5.0*conv, 30.0*conv, 60.0*conv, 30.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 7204083.8568, dp, "Expect d = 7204083.8568 m")
        
        s = boundary_value_method(E, -30.0*conv, 45.0*conv, 30.0*conv, 45.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 6652788.1287, dp, "Expect d = 6652788.1287 m")
        
        s = boundary_value_method(E, -60.0*conv, 60.0*conv, 5.0*conv, 60.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 7213412.4477, dp, "Expect d = 7213412.4477 m")
        
        s = boundary_value_method(E, -75.0*conv, 85.0*conv, 1.0*conv, 85.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 8442938.5899, dp, "Expect d = 8442938.5899 m")
        
        s = boundary_value_method(E, -80.0*conv, 89.5*conv, 0.0, 89.5*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 8888783.7815, dp, "Expect d = 8888783.7815 m")
        
        print("Test passed")
        print(" ")
    
    
    def test_boundary_value_method_table_4(self):
        print(" ")
        print("Test use of the boundary value method on a spheroid, comparing to Table 4 of Panou, 2013")
        a = 6378137.0
        b = 6378137.0   # GRS80 ellipsoid
        c = 6356752.3141
        E = EllipsoidShape(a, b, c)
        conv = math.pi / 180.0 # Conversion factor from degrees to radians
        dp = 4 # Number of decimal places used in comparison
        
        s = boundary_value_method(E, 0.0, 0.0, 0.0, 90.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 10018754.1714, dp, "Expect d =10018754.1714 m")
        
        s = boundary_value_method(E, -1.0*conv, 0.0, 0.0, 179.5*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 19884417.8083, dp, "Expect d = 19884417.8083 m")
        
        s = boundary_value_method(E, 5.0*conv, 0.0, -80.0*conv, 170.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 11652530.7514, dp, "Expect d = 11652530.7514 m")
        
        s = boundary_value_method(E, 30.0*conv, 0.0, -75.0*conv, 120.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 14057886.8752, dp, "Expect d = 14057886.8752 m")
        
        s = boundary_value_method(E, 60.0*conv, 0.0, -60.0*conv, 40.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 13767414.8267, dp, "Expect d = 13767414.8267 m")
        
        s = boundary_value_method(E, 75.0*conv, 0.0, -30.0*conv, 0.5*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 11661713.4496, dp, "Expect d = 11661713.4496 m")
        
        s = boundary_value_method(E, 80.0*conv, 0.0, -5.0*conv, 120.0*conv, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 11105138.2902, dp, "Expect d = 11105138.2902 m")
        
        s = boundary_value_method(E, 0.0, 0.0, 60.0*conv, 0.0, Jacobi=True, n = 16000)      
        self.assertAlmostEqual(s[0], 6663348.2060, dp, "Expect d = 6663348.2060 m")
        
        print("Test passed")
        print(" ")
       
        
            
if __name__ == '__main__':
    unittest.main()


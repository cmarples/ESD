# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 15:21:05 2022

@author: Callum Marples

Test curve tracing of ellipsoid-sphere intersection
"""

import os
os.chdir("..")

import unittest
import math

from leod.newton_raphson import newton_raphson

print(" ")    

class TestCurveTracing(unittest.TestCase):
    

    def test_newton_raphson(self):
        print(" ")
        print("Test Newton-Raphson method on the quadratic x^2 - x - 1 = 0")
        
        fun = lambda x: x**2 - x - 1.0
        fun_prime = lambda x: 2*x - 1.0
        x = newton_raphson(fun, fun_prime, 1.0)
        self.assertAlmostEqual(x, 0.5*(1+math.sqrt(5.0)), 9, "Expect x = phi = 1.61803398874...")
        
        print("Test passed")
        print(" ")

if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 12:36:25 2022

@author: Callum Marples

This file is used for implementing unit tests of the fast marching code.
"""

import os
os.chdir("..")

import unittest
import math

from leod.ellipsoid_shape import EllipsoidShape
from leod.fmm_vertex import FmmVertex

import leod.fmm_polar_graph as pg
import leod.fmm_fast_marching as fm

class TestFMM(unittest.TestCase):
    
    def test_8_neighbour(self):
        print("Simple tests of 8-neighbour GeoGrid")
        a = 3.0
        b = 2.0
        c = 1.0
        shape = EllipsoidShape(a, b, c)
        no_theta = 5
        no_phi = 8
        vertex = pg.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
        
        self.assertEqual(len(vertex), 26, "Should be 26")
        self.assertEqual(vertex[9].neighbour, [1, 17, 16, 10, 8, 2, 24, 18], "Should be [1, 17, 16, 10, 8, 2, 24, 18]")
        self.assertEqual(vertex[8].neighbour, [0, 16, 7, 1, 15, 9], "Should be [0, 16, 7, 1, 15, 9]")
        self.assertEqual(vertex[20].neighbour, [12, 25, 19, 21, 11, 13], "Should be [12, 25, 19, 21, 11, 13]")
        print("Test passed")
        print(" ")
        
    def test_quadratic_terms(self):
        print("Compare alpha, beta and gamma for the two different first order update quadratics")
        shape = EllipsoidShape(1.0, 1.0, 1.0)
        vertex = pg.generate_polar_graph(shape, 19, 36, is_connect_8=True, is_Dijkstra=False)
        fmm = fm.FmmResult(len(vertex))
        self.assertEqual(fmm.accepted[218], False, "Should be False")
        fmm.distance[254] = 0.2455756079379457     # Trial
        fmm.distance[253] = 0.17431148549531636    # Support
        # First order FMM update using the 'trigonometric form'
        [t, alpha, beta, gamma] = fm.fmm_first_order_update_original(218, 254, 253, vertex, fmm)
        [t_lin, alpha_lin, beta_lin, gamma_lin] = fm.fmm_first_order_update(218, 254, 253, vertex, fmm)
        
        beta_modified = 2*alpha_lin*fmm.distance[253] + beta_lin
        gamma_modified = alpha_lin*fmm.distance[253]*fmm.distance[253] + beta_lin*fmm.distance[253] + gamma_lin
        
        
        self.assertAlmostEqual(alpha, alpha_lin, 9, "Should be 0.029468289375376946")
        self.assertAlmostEqual(beta, beta_modified, 9, "Should be -0.004103858786155879")
        self.assertAlmostEqual(gamma, gamma_modified, 9, "Should be -0.0005978125439864344")
        self.assertAlmostEqual(t, t_lin, 9, "Should be 0.40248418407215864")

if __name__ == '__main__':
    unittest.main()
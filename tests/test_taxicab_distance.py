# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 23:19:28 2022

@author: Cal
"""

import os
os.chdir("..")

import unittest

from esd.geo.taxicab import sphere_tcd, spheroid_tcd, triaxial_tcd
from esd.shape import EllipsoidShape
from esd.fmm.mesh_pol import gen_pol_mesh
from esd.fmm.callers import distance_pair

class TestTaxicab(unittest.TestCase):

    def test_taxicab_sphere(self):
        print("Test taxicab distance and 4-neighbour Dijkstra on a sphere")
        start_point = [50.0, 60.0]
        end_point = [90.0, 0.0]
        d = sphere_tcd(1.0, start_point, end_point, is_radians=False)
        # Compare with use of analytical formula
        self.assertAlmostEqual(d, 1.500331566, 9, "Expect d = 1.500331566")
        # Compare with Dijkstra's algorithm on a 4 neighbour polar mesh
        shape = EllipsoidShape(1.0, 1.0, 1.0)
        mesh1 = gen_pol_mesh(91, 180, shape, is_connect_8=False)
        d_mesh, fmm = distance_pair(shape, mesh1, start_point, end_point, is_dijkstra=True, is_radians=False)
        self.assertAlmostEqual(d, d_mesh, 2, "Expect True")
        print("Test passed")
        print(" ")
        
    def test_taxicab_oblate(self):
        print("Test taxicab distance and 4-neighbour Dijkstra on an oblate spheroid")
        start_point = [50.0, 60.0]
        end_point = [90.0, 0.0]
        shape = EllipsoidShape(2.0, 2.0, 1.0)
        d = spheroid_tcd(shape.a_axis, shape.c_axis, start_point, end_point, is_radians=False)
        # Compare with Dijkstra's algorithm on a 4 neighbour polar mesh
        mesh1 = gen_pol_mesh(91, 180, shape, is_connect_8=False)
        d_mesh, fmm1 = distance_pair(shape, mesh1, start_point, end_point, is_dijkstra=True, is_radians=False)
        self.assertAlmostEqual(d, d_mesh, 2, "Expect True")
        
    def test_taxicab_prolate(self):
        print("Test taxicab distance and 4-neighbour Dijkstra on a prolate spheroid")
        start_point = [50.0, 60.0]
        end_point = [90.0, 0.0]
        shape = EllipsoidShape(1.0, 1.0, 2.0)
        d = spheroid_tcd(shape.a_axis, shape.c_axis, start_point, end_point, is_radians=False)
        # Compare with Dijkstra's algorithm on a 4 neighbour polar mesh
        mesh1 = gen_pol_mesh(91, 180, shape, is_connect_8=False)
        d_mesh, fmm = distance_pair(shape, mesh1, start_point, end_point, is_dijkstra=True, is_radians=False)
        self.assertAlmostEqual(d, d_mesh, 2, "Expect True")
        
    def test_taxicab_triaxial(self):
        print("Test taxicab distance and 4-neighbour Dijkstra on a triaxial ellipsoid")
        start_point = [50.0, 60.0]
        end_point = [90.0, 0.0]
        shape = EllipsoidShape(3.0, 2.0, 1.0)
        d = triaxial_tcd(shape, start_point, end_point, is_radians=False)
        # Compare with Dijkstra's algorithm on a 4 neighbour polar mesh
        mesh1 = gen_pol_mesh(91, 180, shape, is_connect_8=False)
        d_mesh, fmm = distance_pair(shape, mesh1, start_point, end_point, is_dijkstra=True, is_radians=False)
        self.assertAlmostEqual(d, d_mesh, 2, "Expect True")

if __name__ == '__main__':
    unittest.main()
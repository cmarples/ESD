# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 23:19:28 2022

@author: Cal
"""

import os
os.chdir("..")

import unittest
from math import sqrt
from esd.shape import EllipsoidShape
from esd.fmm.mesh_pol import gen_pol_mesh
from esd.fmm.mesh_ico import gen_ico_mesh
from esd.fmm.callers import distance_pair, distance_multiple, distance_all, distance_end
from esd.geo.triaxial import bvm_dist
from esd.geo.sphere import gc_dist

class TestFMM(unittest.TestCase):

    def test_fast_marching(self):
        print("Test the fast marching routines")
        ### Input parameters
        # Set the parameters of each grid so their resolutions are as close as possible. 
        no_vertices = 4000
        # Solve quadratic 10*n_d^2 + 2 - n_v = 0 for n_d, the number of face
        # divisions in the icosahedral triangulation.
        no_divisions = round(sqrt(0.1*(no_vertices - 2)))
        # Solve quadratic 0.5*n_phi^2 - n_phi + 2 - n_v = 0 for n_phi, the number
        # of phi values in the theta-phi grid.
        no_phi = round(1 + sqrt(2*no_vertices - 3))
        if no_phi % 2 == 1:
            # Make no_phi even by adding one.
            no_phi += 1
        no_theta = round(0.5*no_phi) + 1
        
        start_point = [50.0, 60.0]
        end_point = [90.0, 0.0]
        
        ### Unit sphere
        shape = EllipsoidShape(1.0, 1.0, 1.0)
        # Fast marching method
        mesh1 = gen_pol_mesh(no_theta, no_phi, shape, is_connect_8=True)
        mesh2 = gen_ico_mesh(no_divisions, shape, is_split=False)
        mesh3 = gen_ico_mesh(no_divisions, shape, is_split=True)
        d1, fmm1 = distance_pair(shape, mesh1, start_point, end_point, is_dijkstra=False, is_radians=False)
        d2, fmm2 = distance_pair(shape, mesh2, start_point, end_point, is_dijkstra=False, is_radians=False)
        d3, fmm3 = distance_pair(shape, mesh3, start_point, end_point, is_dijkstra=False, is_radians=False)
        # For approximate equality, require agreement to one decimal place
        self.assertAlmostEqual(d1, d2, 1, "Expect True")
        self.assertAlmostEqual(d1, d3, 1, "Expect True")
        self.assertAlmostEqual(d2, d3, 1, "Expect True")
        # True distance
        s = gc_dist(1.0, start_point, end_point, is_radians=False)
        self.assertAlmostEqual(d1, s, 1, "Expect True")
        self.assertAlmostEqual(d2, s, 1, "Expect True")
        self.assertAlmostEqual(d3, s, 1, "Expect True")
        
        ### Triaxial ellipsoid
        shape = EllipsoidShape(3.0, 2.0, 1.0)
        # Fast marching method
        mesh1 = gen_pol_mesh(no_theta, no_phi, shape, is_connect_8=True)
        mesh2 = gen_ico_mesh(no_divisions, shape, is_split=False)
        mesh3 = gen_ico_mesh(no_divisions, shape, is_split=True)
        d1, fmm1 = distance_pair(shape, mesh1, start_point, end_point, is_dijkstra=False, is_radians=False)
        d2, fmm2 = distance_pair(shape, mesh2, start_point, end_point, is_dijkstra=False, is_radians=False)
        d3, fmm3 = distance_pair(shape, mesh3, start_point, end_point, is_dijkstra=False, is_radians=False)
        self.assertAlmostEqual(d1, d2, 1, "Expect True")
        self.assertAlmostEqual(d1, d3, 1, "Expect True")
        self.assertAlmostEqual(d2, d3, 1, "Expect True")
        # Boundary value method
        s, s_path = bvm_dist(shape, start_point, end_point, is_radians=False, is_jacobi=False, n=10000)
        self.assertAlmostEqual(d1, s, 1, "Expect True")
        self.assertAlmostEqual(d2, s, 1, "Expect True")
        self.assertAlmostEqual(d3, s, 1, "Expect True")
        
        ### Multiple distance
        ends = [end_point, [120.0, 40.0]]
        dm, fmm4 = distance_multiple(shape, mesh1, start_point, ends, is_dijkstra=False, is_radians=False)
        self.assertEqual(dm[0], d1)
        dm, fmm4 = distance_multiple(shape, mesh2, start_point, ends, is_dijkstra=False, is_radians=False)
        self.assertEqual(dm[0], d2)
        dm, fmm4 = distance_multiple(shape, mesh3, start_point, ends, is_dijkstra=False, is_radians=False)
        self.assertEqual(dm[0], d3)
        
        ### All distance
        fmm_all_1 = distance_all(shape, mesh1, start_point, is_dijkstra=False, is_radians=False)
        d1_all = distance_end(shape, mesh1, fmm_all_1, ends[0], is_radians=False)
        self.assertEqual(d1_all, d1)
        fmm_all_2 = distance_all(shape, mesh2, start_point, is_dijkstra=False, is_radians=False)
        d2_all = distance_end(shape, mesh2, fmm_all_2, ends[0], is_radians=False)
        self.assertEqual(d2_all, d2)
        fmm_all_3 = distance_all(shape, mesh3, start_point, is_dijkstra=False, is_radians=False)
        d3_all = distance_end(shape, mesh3, fmm_all_3, ends[0], is_radians=False)
        self.assertEqual(d3_all, d3)
        
        print("Test passed")
        print(" ")
        
    
if __name__ == '__main__':
    unittest.main()

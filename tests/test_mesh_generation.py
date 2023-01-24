# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:03:10 2023

@author: Cal
"""

import unittest

from math import pi
from esd.shape import EllipsoidShape
from esd.grid import Grid, binary_search
from esd.fmm.mesh_pol import gen_pol_mesh
from esd.fmm.mesh_ico import gen_ico_mesh

class TestBinarySearch(unittest.TestCase):
    
    def test_binary_search(self):
        # Test the binary_search routine by using a simple list example, with
        # a set of points.
        print("Test binary search")
        v_list = [40.0, 41.0, 42.0, 43.0, 44.0, 45.0]
        x1 = 41.8
        v1 = binary_search(v_list, x1)
        self.assertEqual(v1, 2, "Expect 2")
        x2 = 42.4
        v2 = binary_search(v_list, x2)
        self.assertEqual(v2, 2, "Expect 2")
        x3 = 40.1
        v3 = binary_search(v_list, x3)
        self.assertEqual(v3, 0, "Expect 0")
        x4 = 38.9
        v4 = binary_search(v_list, x4)
        self.assertEqual(v4, 0, "Expect 0")
        x5 = 44.9
        v5 = binary_search(v_list, x5)
        self.assertEqual(v5, 5, "Expect 5")
        x6 = 47.0
        v6 = binary_search(v_list, x6)
        self.assertEqual(v6, 5, "Expect 5")
        print("Test passed")
        print(" ")
        
class TestPolar(unittest.TestCase):

    def test_polar_grid(self):
        # Test the polar grid parameters, the functions for converting
        # scalar and vector indices and the functions for finding closest vertices.
        print("Test grid generation")
        grid = Grid(19, 36)
        # Grid parameters.
        self.assertEqual(grid.no_vertices, 614, "Expect 614")
        self.assertEqual(len(grid.theta_list), 19, "Expect 19")
        self.assertEqual(len(grid.phi_list), 37, "Expect 36+1=37")
        # Scalar/vector index conversions.
        self.assertEqual(grid.get_vertex_index(12, 19), 416, "Expect 416")
        self.assertEqual(grid.get_theta_index(416), 12, "Expect 12")
        self.assertEqual(grid.get_phi_index(416, 12), 19, "Expect 19")
        # Finding closest vertices
        th = 116.0 * pi/180.0
        ph = 191.5 * pi/180.0
        self.assertEqual(grid.find_theta_index(th), 12, "Expect 12")
        self.assertEqual(grid.find_phi_index(ph), 19, "Expect 19")
        self.assertEqual(grid.find_vertex_index(th, ph), 416, "Expect 416")
        th = 2.0 * pi/180.0
        ph = 50.0 * pi/180.0
        self.assertEqual(grid.find_theta_index(th), 0, "Expect 0")
        self.assertEqual(grid.find_phi_index(ph), 5, "Expect 5")
        self.assertEqual(grid.find_vertex_index(th, ph), 0, "Expect 0")
        th = 174.9 * pi/180.0
        ph = 359.5 * pi/180.0
        self.assertEqual(grid.find_theta_index(th), 17, "Expect 17")
        self.assertEqual(grid.find_phi_index(ph), 0, "Expect 0")
        self.assertEqual(grid.find_vertex_index(th, ph), 577, "Expect 577")
        print("Test passed")
        print(" ")
        
    def test_polar_mesh(self):
        # Test the 4 and 8 neighbour polar meshes.
        print("Test polar mesh")
        shape = EllipsoidShape(3.0, 2.0, 1.0)
        # 4 neighbour mesh.
        mesh_4 = gen_pol_mesh(9, 8, shape, is_connect_8=False)
        self.assertEqual(mesh_4.no_vertices, 58, "Expect 58")
        # Number of neighbours.
        self.assertEqual(len(mesh_4.vertex[0].neighbour), 8, "Expect 8")
        self.assertEqual(len(mesh_4.vertex[57].neighbour), 8, "Expect 8")
        self.assertEqual(len(mesh_4.vertex[28].neighbour), 4, "Expect 4")
        # Neighbours of given vertex (test phi boundary and lack of diagonals).
        self.assertEqual(25 in mesh_4.vertex[32].neighbour.keys(), True, "Expect True")
        self.assertEqual(33 in mesh_4.vertex[32].neighbour.keys(), False, "Expect False")
  
        # 8 neighbour mesh.
        mesh_8 = gen_pol_mesh(9, 8, shape, is_connect_8=True)
        self.assertEqual(mesh_8.no_vertices, 58, "Expect 58")
        # Number of neighbours.
        self.assertEqual(len(mesh_8.vertex[0].neighbour), 8, "Expect 8")
        self.assertEqual(len(mesh_8.vertex[57].neighbour), 8, "Expect 8")
        self.assertEqual(len(mesh_8.vertex[28].neighbour), 8, "Expect 8")
        # Neighbours of given vertex (test phi boundary and lack of diagonals).
        self.assertEqual(25 in mesh_8.vertex[32].neighbour.keys(), True, "Expect True")
        self.assertEqual(33 in mesh_8.vertex[32].neighbour.keys(), True, "Expect True")
        
        # Neighbour distances
        self.assertAlmostEqual(mesh_4.vertex[0].neighbour[1].distance, 1.15057108, 8, "Expect 1.15057108")
        mesh_1 = gen_pol_mesh(81, 80, shape, is_connect_8=True)
        for neigh in mesh_1.vertex[0].neighbour.values():
            self.assertLess(neigh.distance, 0.2, "Expect small distance")
        for neigh in mesh_1.vertex[1].neighbour.values():
            self.assertLess(neigh.distance, 0.2, "Expect small distance")
        for neigh in mesh_1.vertex[250].neighbour.values():
            self.assertLess(neigh.distance, 0.2, "Expect small distance")
        self.assertEqual(mesh_1.vertex[0].neighbour[1].distance, mesh_1.vertex[1].neighbour[0].distance)
        
        print("Test passed")
        print(" ")

class TestIco(unittest.TestCase):
    
    def test_ico_mesh(self):
        # Test the icosahedral mesh.
        print("Test icosahedral mesh")
        shape = EllipsoidShape(3.0, 2.0, 1.0)
        mesh_ico = gen_ico_mesh(4, shape, is_split=False)
        # Number of vertices.
        self.assertEqual(len(mesh_ico.vertex), 162, "Expect 162")
        # Number of neighbours.
        self.assertEqual(len(mesh_ico.vertex[0].neighbour), 5, "Expect 5")
        self.assertEqual(len(mesh_ico.vertex[1].neighbour), 6, "Expect 6")
        # Number of faces and edges.
        no_faces = 0
        no_edges = 0
        for i in range(len(mesh_ico.vertex)):
            no_faces += len(mesh_ico.vertex[i].face) # This will triple count faces.
            no_edges += len(mesh_ico.vertex[i].neighbour.keys()) # This will double count edges.
        self.assertEqual(no_faces/3, 320, "Expect 320")
        self.assertEqual(no_edges/2, 480, "Expect 480")
        
        # Neighbour distances
        mesh_1 = gen_ico_mesh(15, shape, is_split=False)
        for neigh in mesh_1.vertex[0].neighbour.values():
            self.assertLess(neigh.distance, 0.2, "Expect small distance")
        for neigh in mesh_1.vertex[1].neighbour.values():
            self.assertLess(neigh.distance, 0.2, "Expect small distance")
        for neigh in mesh_1.vertex[250].neighbour.values():
            self.assertLess(neigh.distance, 0.2, "Expect small distance")
        self.assertEqual(mesh_1.vertex[0].neighbour[1].distance, mesh_1.vertex[1].neighbour[0].distance)
        
        print("Test passed")
        print(" ")
      
if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:54:00 2022

@author: Callum Marples

This file is used for devopment purposes (to test/debug the code).
It works by generating the necessary objects and calculating the distances
for a simple example.
"""

import os
os.chdir("..")

import math
import time

from leod.ellipsoid_shape import EllipsoidShape
from leod.fmm_vertex import FmmVertex
from leod.fmm_polar_graph import generate_polar_graph
from leod.fmm_polar_graph import check_triangles
from leod.fmm_polar_graph import find_vertex_index
from leod.fmm_fast_marching import fast_marching
from leod.sphere_geodesics import great_circle_distance

test_no = 2

if test_no == 1: # Create ellipsoid grid using polar coordinates
    shape = EllipsoidShape(1.0, 1.0, 1.0)
    no_theta = 181
    no_phi = 360
    vertex = generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    
    # Check that these triangles are acute (and thus valid for the FMM)
    #no_obtuse = check_triangles(vertex)
    
    # Structured grid information
    no_vertices = (no_theta - 2)*no_phi + 2
    delta_theta = math.pi / (no_theta - 1)
    delta_phi = 2.0*math.pi / no_phi
    
    # Compute lists of theta and phi values
    theta_list = [0.0] * no_theta
    phi_list = [0.0] * (no_phi+1)
    for i in range(no_theta):
        theta_list[i] = i *delta_theta
    for i in range(no_phi+1):
        phi_list[i] = i * delta_phi
    
    start_th = 90.0 * math.pi / 180.0
    start_ph = 0.0  * math.pi / 180.0  
    end_th = 50.0 * math.pi / 180.0
    end_ph = 60.0 * math.pi / 180.0
    start_vertex = find_vertex_index(theta_list, phi_list, start_th, start_ph)
    end_vertex = find_vertex_index(theta_list, phi_list, end_th, end_ph)
    
    fmm = fast_marching(vertex, start_vertex, 1)
    d = fmm.distance[end_vertex]
    
    s = great_circle_distance(1.0, start_th, start_ph, end_th, end_ph)
    
if test_no == 2: # Test triaxial graph for obtuse angles
    shape = EllipsoidShape(3.0, 2.0, 1.0)
    no_theta = 181
    no_phi = 360
    vertex = generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    
    # Check that these triangles are acute (and thus valid for the FMM)
    no_obtuse = check_triangles(vertex)
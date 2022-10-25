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

import leod.fmm_polar_graph as pg
from leod.ellipsoid_shape import EllipsoidShape
from leod.fmm_vertex import FmmVertex
from leod.fmm_polar_graph import generate_polar_graph
from leod.fmm_polar_graph import check_triangles
from leod.fmm_polar_graph import find_vertex_index
from leod.fmm_fast_marching import fast_marching
from leod.fmm_fast_marching import fmm_endpoint
from leod.fmm_fast_marching import fmm_idw
from leod.sphere_geodesics import great_circle_distance
from leod.triaxial_geodesics import boundary_value_method

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
    shape = EllipsoidShape(1.0, 1.0, 1.0)
    shape.normalise()
    no_theta = 181
    no_phi = 360
    vertex = generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    
    # Check that these triangles are acute (and thus valid for the FMM)
    no_obtuse_1 = check_triangles(vertex)
    p = pg.find_obtuse_angles(vertex, no_theta, no_phi, len(vertex))
    
    if no_obtuse_1 > 0:
        pg.split_update_triangles(vertex, no_theta, no_phi, len(vertex))
    
    no_obtuse_2 = check_triangles(vertex)
    
    no_removals = no_obtuse_1 - no_obtuse_2
    
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
    
    deg2rad = math.pi / 180.0
    start_th = 90.0 * deg2rad
    start_ph = 0.0  * deg2rad  
    end_th = 50.0 * deg2rad
    end_ph = 60.0 * deg2rad
    start_vertex = find_vertex_index(theta_list, phi_list, start_th, start_ph)
    end_vertex = find_vertex_index(theta_list, phi_list, end_th, end_ph)
    
    order = 1
    fmm = fast_marching(vertex, start_vertex, order)
    
    end_carts = shape.polar2cart(end_th, end_ph)
    if vertex[end_vertex].carts[0] != end_carts[0] or vertex[end_vertex].carts[1] != end_carts[1] or vertex[end_vertex].carts[2] != end_carts[2]:
        temp = fmm.update[end_vertex]
        if isinstance(temp, int) == True:
            trial = temp
            support = -1
        else:
            [trial, support] = temp
        d = fmm_endpoint(end_vertex, end_carts, vertex, fmm, order, end_vertex, trial, support=-1)
        
        # IDW Interpolation
        d2 = fmm_idw(end_carts, vertex[17701].carts, vertex[17341].carts, vertex[17342].carts, vertex[17702].carts,
                     fmm.distance[17701], fmm.distance[17341], fmm.distance[17342], fmm.distance[17702])
        
    else:
        d = fmm.distance[end_vertex]
    
    if shape.a_axis == shape.b_axis and shape.b_axis == shape.c_axis:
        s = great_circle_distance(shape.a_axis, start_th, start_ph, end_th, end_ph)
    else:
        s = boundary_value_method(shape, start_th, start_ph, end_th, end_ph, tol=1e-12, Jacobi=False, n = 20000)
    
if test_no == 3: # Test triaxial graph for obtuse angles - at the pole
    shape = EllipsoidShape(3.0, 2.0, 1.0)
    shape.normalise()
    no_theta = 181
    no_phi = 360
    vertex = generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    
    # Check that these triangles are acute (and thus valid for the FMM)
    no_obtuse_1 = check_triangles(vertex)
    p = pg.find_obtuse_angles(vertex, no_theta, no_phi, len(vertex))
    
    if no_obtuse_1 > 0:
        pg.split_update_triangles(vertex, no_theta, no_phi, len(vertex))
    
    no_obtuse_2 = check_triangles(vertex)
    
    no_removals = no_obtuse_1 - no_obtuse_2
    
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
    
    deg2rad = math.pi / 180.0
    start_th = 1.0 * deg2rad
    start_ph = 0.0  * deg2rad  
    end_th = 5.0 * deg2rad
    end_ph = 180.0 * deg2rad
    start_vertex = find_vertex_index(theta_list, phi_list, start_th, start_ph)
    end_vertex = find_vertex_index(theta_list, phi_list, end_th, end_ph)
    
    fmm = fast_marching(vertex, start_vertex, 1)
    d = fmm.distance[end_vertex]
    
    #s = boundary_value_method(shape, 5.0*deg2rad, 0.0*deg2rad, 5.0*deg2rad, 179.0*deg2rad, tol=1e-12, Jacobi=False, n = 20000)
    
    vertex2 = generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    fmm2 = fast_marching(vertex2, start_vertex, 1)
    d2 = fmm2.distance[end_vertex]
    
    x = fmm.distance[1:361]
    y = fmm2.distance[1:361]
    for i in range(360):
        z = x[i] - y[i]
        if z != 0.0:
            print(i)
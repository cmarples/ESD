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
import numpy as np

from leod.ellipsoid_shape import EllipsoidShape
from leod.fmm_vertex import FmmVertex
from leod.fmm_fast_marching import fast_marching
from leod.fmm_fast_marching import fmm_idw
from leod.fmm_fast_marching import fmm_idw3
from leod.fmm_fast_marching import endpoint_distance
from leod.sphere_geodesics import great_circle_distance
from leod.triaxial_geodesics import boundary_value_method
from leod.triangulation_sphere import triangulate_sphere

from leod.fmm_precalculation import precalculate_grid
from leod.fmm_precalculation import split_obtuse_angles
from leod.fmm_precalculation import find_obtuse_angles

import leod.triangulation as tri
import leod.fmm_polar_graph as pg

from leod.fmm_polar_graph import PolarGrid

from leod.fmm_callers import calculate_pair_distance

test_no = 5

if test_no == 1: # FMM on sphere triangulation
    
    r = 1.0 # radius
    shape = EllipsoidShape(1.0, 1.0, 1.0)
    n = 100  # number of triangular divisions
    tic = time.perf_counter()
    vertex = triangulate_sphere(r, n)
    toc = time.perf_counter()
    print(toc - tic)
    no_obtuse = tri.check_triangles(vertex)
    
    # Closest point to given start
    start_th = 90.0 * math.pi / 180.0
    start_ph = 0.0  * math.pi / 180.0  
    p_start = shape.polar2cart(start_th, start_ph)
    p_start = np.array(p_start)
    end_th = 50.0 * math.pi / 180.0
    end_ph = 60.0 * math.pi / 180.0
    p_end = shape.polar2cart(end_th, end_ph)
    p_end = np.array(p_end)
    
    tic = time.perf_counter()
    v_start = tri.find_closest_vertex(vertex, p_start)
    v_end = tri.find_closest_vertex(vertex, p_end)
    toc = time.perf_counter()
    print(toc - tic)
    
    tic = time.perf_counter()
    fmm = fast_marching(vertex, v_start, p_start, 1)
    toc = time.perf_counter()
    print(toc - tic)
    
    tic = time.perf_counter()
    end_face = tri.find_face(vertex, p_end, v_end)
    d = fmm_idw3( p_end, vertex[end_face[0]].carts, vertex[end_face[1]].carts, vertex[end_face[2]].carts,
                  fmm.distance[end_face[0]], fmm.distance[end_face[1]], fmm.distance[end_face[2]] )
    #carts = [vertex[v_end].carts]
    #dist = [fmm.distance[v_end]]
    #for i in vertex[v_end].distance_to_neighbour.keys():
    #    dist.append(fmm.distance[i])
    #    carts.append(vertex[i].carts)
    #d = fmm_idw(p_end, carts, dist)
    toc = time.perf_counter()
    print(toc - tic)
    
    s = great_circle_distance(1.0, start_th, start_ph, end_th, end_ph)
    
elif test_no == 2: # Test anisotropic scaling of the icosahedral triangulation
    r = 1.0 # radius
    shape = EllipsoidShape(3.0, 2.0, 1.0)
    n = 50  # number of triangular divisions
    tic = time.perf_counter()
    vertex = triangulate_sphere(r, n)
    toc = time.perf_counter()
    print(toc - tic)
    
    # Scale vertices
    for i in range(len(vertex)):
        vertex[i].carts[0] *= shape.a_axis
        vertex[i].carts[1] *= shape.b_axis
        vertex[i].carts[2] *= shape.c_axis
        
    no_obtuse = pg.check_triangles(vertex)
    
elif test_no == 3: # Test triaxial graph for obtuse angles - at the pole
    shape = EllipsoidShape(3.0, 2.0, 1.0)
    shape.normalise()
    no_theta = 181
    no_phi = 360
    vertex = pg.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    
    # Check that these triangles are acute (and thus valid for the FMM)
    no_obtuse_1 = pg.check_triangles(vertex)
    p = pg.find_obtuse_angles(vertex, no_theta, no_phi, len(vertex))
    
    if no_obtuse_1 > 0:
        pg.split_update_triangles(vertex, no_theta, no_phi, len(vertex))
    
    no_obtuse_2 = pg.check_triangles(vertex)
    
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
    start_th = 10.0 * deg2rad
    start_ph = 0.0  * deg2rad  
    end_th = 5.0 * deg2rad
    end_ph = 165.0 * deg2rad
    start_vertex = pg.find_vertex_index(theta_list, phi_list, start_th, start_ph)
    end_vertex = pg.find_vertex_index(theta_list, phi_list, end_th, end_ph)
    start_point = shape.polar2cart(start_th, start_ph)
    
    fmm = fast_marching(vertex, start_vertex, start_point, 1)
    d = fmm.distance[end_vertex]
    
    s = boundary_value_method(shape, start_th, start_ph, end_th, end_ph, tol=1e-12, Jacobi=False, n = 20000)
    
    vertex2 = pg.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    fmm2 = fast_marching(vertex2, start_vertex, start_point, 1)
    d2 = fmm2.distance[end_vertex]
    
    x = fmm.distance[1:361]
    y = fmm2.distance[1:361]
    for i in range(360):
        z = x[i] - y[i]
        if z != 0.0:
            print(i)
            
elif test_no == 4: # Speed of FMM
    shape = EllipsoidShape(3.0, 2.0, 1.0)
    shape.normalise()
    no_theta = 181
    no_phi = 360
    tic = time.perf_counter()
    
    vertex, grid = pg.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    
    toc = time.perf_counter()
    print(toc - tic)
    
    tic = time.perf_counter()
    max_angle = precalculate_grid(vertex)
    toc = time.perf_counter()
    print(toc - tic)

    # Check that these triangles are acute (and thus valid for the FMM)
    #no_obtuse_1 = pg.check_triangles(vertex)
    
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
    start_th = 70.0 * deg2rad
    start_ph = 80.0  * deg2rad  
    end_th = 40.0 * deg2rad
    end_ph = 340.0 * deg2rad
    start_vertex, st_th, st_ph = pg.find_vertex_index(theta_list, phi_list, start_th, start_ph)
    end_vertex, end_th_index, end_ph_index = pg.find_vertex_index(theta_list, phi_list, end_th, end_ph)
    start_point = shape.polar2cart(start_th, start_ph)
    
    end_dict = {}
    end_dict[end_vertex] = False
    for j in vertex[end_vertex].neighbour.keys():
        end_dict[j] = False
    
    tic = time.perf_counter()
    
    fmm = fast_marching(vertex, start_vertex, start_point, 1, end_dict)
    #d = fmm.distance[end_vertex]
    d = endpoint_distance(vertex, fmm, end_th, end_ph, end_vertex, shape)
    
    toc = time.perf_counter()
    print(toc - tic)
    
    s = boundary_value_method(shape, start_th, start_ph, end_th, end_ph, tol=1e-12, Jacobi=False, n = 20000)


elif test_no == 5: # Test pair routine
    shape = EllipsoidShape(1.0, 1.0, 1.0)
    shape.normalise()
    no_theta = 181
    no_phi = 360
    tic = time.perf_counter()
    
    vertex, polar_grid = pg.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=False)
    
    toc = time.perf_counter()
    t1 = toc - tic
    print(t1)
    
    tic = time.perf_counter()
    precalculate_grid(vertex)
    toc = time.perf_counter()
    t2 = toc - tic
    print(t2)

    start_th = 50.0
    start_ph = 60.0  
    end_th = 90.0
    end_ph = 0.0
    
    tic = time.perf_counter()
    
    d, fmm = calculate_pair_distance(shape, vertex, [start_th, start_ph], [end_th, end_ph],
                                     0, graph_type="polar", grid=polar_grid, is_radians=False)
    
    toc = time.perf_counter()
    t3 = toc - tic
    print(t3)
    
    tic = time.perf_counter()
    conv = math.pi / 180.0
    if shape.is_sphere() == True:
        s = great_circle_distance(shape.a_axis, start_th*conv, start_ph*conv, end_th*conv, end_ph*conv)
    else:
        s = boundary_value_method(shape, start_th*conv, start_ph*conv, end_th*conv, end_ph*conv, tol=1e-12, Jacobi=False, n = 20000)
    toc = time.perf_counter()
    t4 = toc - tic
    print(t4)
    
    print('t_FMM = ', t1+t2+t3)
    
elif test_no == 6: # Test pair routine (triangulation)
    shape = EllipsoidShape(3.0, 2.0, 1.0)
    shape.normalise()
    
    tic = time.perf_counter()
    
    n = 100  # number of triangular divisions
    vertex = triangulate_sphere(1.0, n)
    # Scale to ellipsoid
    for i in range(len(vertex)):
        vertex[i].carts[0] *= shape.a_axis
        vertex[i].carts[1] *= shape.b_axis
        vertex[i].carts[2] *= shape.c_axis
    
    toc = time.perf_counter()
    print(toc - tic)
    
    tic = time.perf_counter()
    max_angle = precalculate_grid(vertex)
    [ob1, m1] = find_obtuse_angles(vertex)
    split_obtuse_angles(vertex)
    [ob2, m2] = find_obtuse_angles(vertex)
    toc = time.perf_counter()
    print(toc - tic)

    start_th = 50.0
    start_ph = 60.0  
    end_th = 90.0
    end_ph = 0.0
    
    tic = time.perf_counter()
    
    d, fmm = calculate_pair_distance(shape, vertex, [start_th, start_ph], [end_th, end_ph],
                                     1, graph_type="tri", grid=-1, is_radians=False)
    
    toc = time.perf_counter()
    print(toc - tic)
    
    tic = time.perf_counter()
    conv = math.pi / 180.0
    if shape.is_sphere() == True:
        s = great_circle_distance(shape.a_axis, start_th*conv, start_ph*conv, end_th*conv, end_ph*conv)
    else:
        s = boundary_value_method(shape, start_th*conv, start_ph*conv, end_th*conv, end_ph*conv, tol=1e-12, Jacobi=False, n = 20000)
    toc = time.perf_counter()
    print(toc - tic)
    
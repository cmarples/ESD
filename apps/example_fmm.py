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
from leod.fmm_fast_marching import fmm_endpoint
from leod.fmm_fast_marching import fmm_idw
from leod.sphere_geodesics import great_circle_distance
from leod.triaxial_geodesics import boundary_value_method
from leod.triangulation_sphere import triangulate_sphere

import leod.triangulation as tri

test_no = 1

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
    
    fmm = fast_marching(vertex, v_start, 1)
    d = fmm.distance[v_end]
    
    s = great_circle_distance(1.0, start_th, start_ph, end_th, end_ph)
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
from leod.fmm_fast_marching import fast_marching
from leod.fmm_fast_marching import fmm_endpoint
from leod.fmm_fast_marching import fmm_idw
from leod.sphere_geodesics import great_circle_distance
from leod.triaxial_geodesics import boundary_value_method
from leod.triangulation_sphere import triangulate_sphere
from leod.triangulation import check_triangles

test_no = 1

if test_no == 1: # FMM on sphere triangulation
    
    r = 1.0 # radius
    n = 4  # number of triangular divisions
    tic = time.perf_counter()
    vertex = triangulate_sphere(r, n)
    toc = time.perf_counter()
    print(toc - tic)
    no_obtuse = check_triangles(vertex)
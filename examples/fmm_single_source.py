# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 12:38:50 2022

@author: Callum Marples

Example : single-source to many endpoint problem using the fast marching method.
"""

import os
os.chdir("..")

import math
import time
import numpy as np

import leod

### Define shape
shape = leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 1.0)
shape.normalise()

### Start and end points
start_point = [90.0, 0.0]
end_point = [ [50.0, 60.0], [0.0, 0.0] ]

### Polar mesh
no_theta = 51
no_phi = 100
print('')
print('Run the fast marching method on a polar mesh.')

tic = time.perf_counter()
vertex, polar_grid = leod.fmm_polar_graph.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True)
toc = time.perf_counter()
t1 = toc - tic
print('Run-time (generate_polar_graph) =', t1)

tic = time.perf_counter()
leod.fmm_precalculation.precalculate_grid(vertex)
toc = time.perf_counter()
t2 = toc - tic
print('Run-time (precalculate_grid) =', t2)

### Fast marching on the polar mesh
tic = time.perf_counter()
d, fmm = leod.fmm_callers.calculate_distances(shape, vertex, start_point, end_point, 1,
                                              graph_type="polar", grid=polar_grid, is_radians=False)
toc = time.perf_counter()
t3 = toc - tic
print('Run-time (fast marching) =', t3)


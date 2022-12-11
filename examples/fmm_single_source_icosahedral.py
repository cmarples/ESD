# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 17:43:21 2022

@author: Callum Marples
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
start_point = [0.0, 0.0]
end_point = [180.0, 0.0]

### Icosahedral grid
no_divisions = 77
print('')
print('Run the fast marching method on an icosahedral mesh.')

tic = time.perf_counter()
vertex = leod.triangulation_sphere.triangulate_sphere(1.0, no_divisions)
# Scale to ellipsoid
tic = time.perf_counter()
for k in range(len(vertex)):
    vertex[k].carts[0] *= shape.a_axis
    vertex[k].carts[1] *= shape.b_axis
    vertex[k].carts[2] *= shape.c_axis
toc = time.perf_counter()
t1 = toc - tic
print('Run-time (generate_icosahedral_graph) =', t1)

tic = time.perf_counter()
leod.fmm_precalculation.precalculate_grid(vertex)
#leod.fmm_precalculation.split_obtuse_angles(vertex)
toc = time.perf_counter()
t2 = toc - tic
print('Run-time (precalculate_grid) =', t2)

### Fast marching on the polar mesh
tic = time.perf_counter()
d, fmm = leod.fmm_callers.calculate_pair_distance(shape, vertex, start_point, end_point,
                                                  1, graph_type="tri", grid=-1, is_radians=False)
toc = time.perf_counter()
t3 = toc - tic
print('Run-time (fast marching) =', t3)
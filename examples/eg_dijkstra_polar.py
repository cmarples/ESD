# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 19:18:48 2022

@author: Cal
"""

import leod
import math
import numpy as np
import time
import os


os.chdir("..")

### Define shape
shape = leod.shape.EllipsoidShape(1.0, 1.0, 1.0)
shape.normalise()

### Start and end points
start_point = [50.0, 60.0]
end_point = [90.0, 0.0]

### Polar mesh
no_theta = 181
no_phi = 360
print("")
print("Run the fast marching method on a polar mesh.")

tic = time.perf_counter()
grid = leod.fmm.grid_pol.gen_pol_grid(no_theta, no_phi, shape, is_connect_8=False)
toc = time.perf_counter()
t1 = toc - tic
print("Run-time (generate_polar_graph) =", t1)

### Dijkstra's algorithm on the polar grid
tic = time.perf_counter()
d, fmm = leod.fmm.callers.distance_pair(shape, grid, start_point, end_point,
                                        is_dijkstra=False, is_radians=False)
toc = time.perf_counter()
t2 = toc - tic
print("Run-time (fast marching) =", t2)
print("Distance =", d)
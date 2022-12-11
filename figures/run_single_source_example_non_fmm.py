# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 23:54:08 2022

@author: Callum Marples

Generate systematic single-source distance data using the true distance 
formula (sphere) and the boundary value method (ellipsoid).

This creates data for Figures 8 and 9 of the Geodesic Paper.
"""

import leod
import numpy as np
import math
import time
import csv
import os

save_flag = True

os.chdir("..")

### Start and end points
start_point = [90.0, 0.0]

# Use a polar grid to define a set of end points.
n_th = 19
n_ph = 18
n_ends = (n_th - 2) * n_ph + 2

vals_theta = np.linspace(0.0, 180.0, num=n_th)
vals_phi = np.linspace(0.0, 360.0, num=n_ph+1)
vals_phi = vals_phi[0:len(vals_phi)-1]

end_point = [0] * n_ends
k = 0
for i in range(n_th):
    if vals_theta[i] == 0.0 or vals_theta[i] == 180.0:
        end_point[k] = [vals_theta[i], 0.0]
        k += 1
    else:
        for j in range(n_ph):
            end_point[k] = [vals_theta[i], vals_phi[j]]
            k += 1

conv = math.pi / 180.0 # Degrees to radians.
d_sphere = [0] * n_ends
t_sphere = [0] * n_ends
d_ellipsoid = [0] * n_ends
t_ellipsoid = [0] * n_ends

### Sphere data
shape1 = leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 1.0)
for i in range(n_ends):
    tic = time.perf_counter()
    d_sphere[i] = leod.sphere_geodesics.great_circle_distance(1.0, start_point[0]*conv, start_point[1]*conv,
                                                              end_point[i][0]*conv, end_point[i][1]*conv)
    toc = time.perf_counter()
    t_sphere[i] = toc - tic

### Ellipsoid data
shape2 = leod.ellipsoid_shape.EllipsoidShape(3.0, 2.0, 1.0)
shape2.normalise()
for i in range(n_ends):
    try:
        tic = time.perf_counter()
        d = leod.triaxial_geodesics.boundary_value_method(shape2, start_point[0]*conv, start_point[1]*conv,
                                                          end_point[i][0]*conv, end_point[i][1]*conv,
                                                          tol=1e-12, Jacobi=False, n = 20000)
        d_ellipsoid[i] = d[0]
        toc = time.perf_counter()
        t_ellipsoid[i] = toc - tic
    except:
        d_ellipsoid[i] = -1.0
        t_ellipsoid[i] = -1.0
        
    if math.isnan(d_ellipsoid[i]):
        d_ellipsoid[i] = -1.0

### Write to files
file_name_s = 'data/geodesics/single_source_sphere_true.txt'
with open(file_name_s, mode="w", newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for k in range(len(d_sphere)):
        writer.writerow([end_point[k][0], end_point[k][1], d_sphere[k], t_sphere[k]])

file_name_e = 'data/geodesics/single_source_ellipsoid_bvm.txt'
with open(file_name_e, mode="w", newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for k in range(len(d_ellipsoid)):
        writer.writerow([end_point[k][0], end_point[k][1], d_ellipsoid[k], t_ellipsoid[k]])
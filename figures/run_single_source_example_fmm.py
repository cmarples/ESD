# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 20:28:39 2022

@author: Callum Marples

Generate systematic single-source distance data using the fast marching method
on a polar mesh and an icosahedral triangulation.

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

### Generate meshes

# Shapes
shape1 = leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 1.0)
shape2 = leod.ellipsoid_shape.EllipsoidShape(3.0, 2.0, 1.0)
shape2.normalise()

# Mesh resolution
no_vertices = 60000

# Solve quadratic 10*n_d^2 + 2 - n_v = 0 for n_d, the number of face
# divisions in the icosahedral triangulation.
no_divisions = round(math.sqrt(0.1*(no_vertices - 2)))
no_vertices_icosahedral = 10*no_divisions**2 + 2

# Solve quadratic 0.5*n_phi^2 - n_phi + 2 - n_v = 0 for n_phi, the number
# of phi values in the theta-phi grid.
no_phi = round(1 + math.sqrt(2*no_vertices - 3))
if no_phi % 2 == 1:
    # Make no_phi even by adding one.
    no_phi += 1
no_theta = round(0.5*no_phi) + 1
no_vertices_polar = no_phi*(no_theta-2) + 2

# Polar meshes
tic = time.perf_counter()
vertex_polar1, polar_grid1 = leod.fmm_polar_graph.generate_polar_graph(shape1, no_theta, no_phi, is_connect_8=True)
leod.fmm_precalculation.precalculate_grid(vertex_polar1)
toc = time.perf_counter()
t_polar1_gen = toc - tic
print('Run-time (sphere polar generation) =', t_polar1_gen)

tic = time.perf_counter()
vertex_polar2, polar_grid2 = leod.fmm_polar_graph.generate_polar_graph(shape2, no_theta, no_phi, is_connect_8=True)
leod.fmm_precalculation.precalculate_grid(vertex_polar2)
toc = time.perf_counter()
t_polar2_gen = toc - tic
print('Run-time (triaxial polar generation) =', t_polar2_gen)

# Icosahedral meshes
tic = time.perf_counter()
vertex_ico1 = leod.triangulation_sphere.triangulate_sphere(1.0, no_divisions)
leod.fmm_precalculation.precalculate_grid(vertex_ico1)
toc = time.perf_counter()
t_ico1_gen = toc - tic
print('Run-time (sphere icosahedral triangulation generation) =', t_ico1_gen)

tic = time.perf_counter()
vertex_ico2 = leod.triangulation_sphere.triangulate_sphere(1.0, no_divisions)
for k in range(len(vertex_ico2)):
    vertex_ico2[k].carts[0] *= shape2.a_axis
    vertex_ico2[k].carts[1] *= shape2.b_axis
    vertex_ico2[k].carts[2] *= shape2.c_axis
leod.fmm_precalculation.precalculate_grid(vertex_ico2)
leod.fmm_precalculation.split_obtuse_angles(vertex_ico2)                            
toc = time.perf_counter()
t_ico2_gen = toc - tic
print('Run-time (triaxial icosahedral triangulation generation) =', t_ico2_gen)
                            
### Sphere data (polar)
tic = time.perf_counter()
d_polar1, fmm_polar1 = leod.fmm_callers.calculate_distances(shape1, vertex_polar1, start_point, end_point, 1,
                                                            graph_type="polar", grid=polar_grid1)
toc = time.perf_counter()
t_polar1_fmm = toc - tic
print('Run-time (sphere polar fmm) =', t_polar1_fmm)

### Sphere data (icosahedral)
tic = time.perf_counter()
d_ico1, fmm_ico1 = leod.fmm_callers.calculate_distances(shape1, vertex_ico1, start_point, end_point, 1)
toc = time.perf_counter()
t_ico1_fmm = toc - tic
print('Run-time (sphere icosahedral triangulation fmm) =', t_ico1_fmm)

### Ellipsoid data (polar)
tic = time.perf_counter()
d_polar2, fmm_polar2 = leod.fmm_callers.calculate_distances(shape2, vertex_polar2, start_point, end_point, 1,
                                                            graph_type="polar", grid=polar_grid2)
toc = time.perf_counter()
t_polar2_fmm = toc - tic
print('Run-time (triaxial polar fmm) =', t_polar2_fmm)

### Ellipsoid data (icosahedral)
tic = time.perf_counter()
d_ico2, fmm_ico2 = leod.fmm_callers.calculate_distances(shape2, vertex_ico2, start_point, end_point, 1)
toc = time.perf_counter()
t_ico2_fmm = toc - tic
print('Run-time (triaxial icosahedral triangulation fmm) =', t_ico2_fmm)

### Write to file
file_name = 'data/geodesics/single_source_fast_marching.txt'
with open(file_name, mode="w", newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for k in range(n_ends):
        writer.writerow([d_polar1[k], d_ico1[k], d_polar2[k], d_ico2[k]])
    writer.writerow([t_polar1_gen, t_ico1_gen, t_polar2_gen, t_ico2_gen])
    writer.writerow([t_polar1_fmm, t_ico1_fmm, t_polar2_fmm, t_ico2_fmm])

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:18:06 2022

@author: Callum Marples

Calculate distance for a particular example, at a range of resolutions.

Example: (theta, phi) = (50, 60) -> (90, 0)
(above angles in degrees)

Five different shapes are used and for each, an appropriate alternate method
is used for comparison.
"""

import leod
import math
import numpy as np
import time
import copy
import csv
import sys
import os

modulename = 'leod.spheroid_geodesics'
if modulename in sys.modules:
    geo_flag = True
else:
    geo_flag = False

os.chdir("..")

### Define shapes.
shape = [ leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 1.0),   # Sphere.
          leod.ellipsoid_shape.EllipsoidShape(6378137.0, 6378137.0, 6356752.3142),  # Earth (WGS84 ellipsoid).
          leod.ellipsoid_shape.EllipsoidShape(2.0, 2.0, 1.0),   # Oblate.
          leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 2.0),   # Prolate.
          leod.ellipsoid_shape.EllipsoidShape(3.0, 2.0, 1.0) ]  # Triaxial.
shape_name = ['sphere', 'wgs84', 'oblate', 'prolate', 'triaxial']

# Normalise shapes so that (abc) = 1.
# Unit sphere already normalised, do not want to do this for the WGS84 example.
shape[2].normalise()
shape[3].normalise()
shape[4].normalise()

# Define start and end points of the geodesic example.
start_point = [50.0, 60.0]
end_point = [90.0, 0.0]
conv = math.pi / 180.0 # Degrees to radians.

# Number of vertices (logarithmic scale)
log_n = np.linspace(2, 5, 25)
no_vertices_float = np.power(10.0, log_n)
no_vertices = [0] * len(no_vertices_float)
for i in range(len(no_vertices)):
    no_vertices[i] = int(round(no_vertices_float[i]))

### Find mesh parameters from the number of vertices.
no_vertices_polar = [0] * len(no_vertices)
no_vertices_icosahedral = [0] * len(no_vertices)
no_divisions = [0] * len(no_vertices)
no_phi = [0] * len(no_vertices)
no_theta = [0] * len(no_vertices)
for i in range(len(no_vertices)):
    # Solve quadratic 10*n_d^2 + 2 - n_v = 0 for n_d, the number of face
    # divisions in the icosahedral triangulation.
    no_divisions[i] = round(math.sqrt(0.1*(no_vertices[i] - 2)))
    no_vertices_icosahedral[i] = 10*no_divisions[i]**2 + 2
    # Solve quadratic 0.5*n_phi^2 - n_phi + 2 - n_v = 0 for n_phi, the number
    # of phi values in the theta-phi grid.
    no_phi[i] = round(1 + math.sqrt(2*no_vertices[i] - 3))
    if no_phi[i] % 2 == 1:
        # Make no_phi even by adding one.
        no_phi[i] += 1
    no_theta[i] = round(0.5*no_phi[i]) + 1
    no_vertices_polar[i] = no_phi[i]*(no_theta[i]-2) + 2

# Write number of vertices data
file_name = 'data/resolution/resolution_number_of_vertices.txt'
with open(file_name, mode="w", newline='') as f:
    
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for i in range(len(no_vertices)):
        row = [no_vertices[i], no_vertices_polar[i], no_vertices_icosahedral[i]]
        writer.writerow(row)
        
# Initialise result lists
t_polar_gen_4 = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_4 = np.zeros([len(no_vertices), len(shape)])
t_fmm_4 = np.zeros([len(no_vertices), len(shape)])
d_dijkstra_4 = np.zeros([len(no_vertices), len(shape)])
d_fmm_4 = np.zeros([len(no_vertices), len(shape)])

t_polar_gen_8 = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_8 = np.zeros([len(no_vertices), len(shape)])
t_fmm_8 = np.zeros([len(no_vertices), len(shape)])
d_dijkstra_8 = np.zeros([len(no_vertices), len(shape)])
d_fmm_8 = np.zeros([len(no_vertices), len(shape)])

t_tri_gen = [0.0] * len(no_vertices)

t_tri_scale_gen = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_tri = np.zeros([len(no_vertices), len(shape)])
t_fmm_tri = np.zeros([len(no_vertices), len(shape)])
d_dijkstra_tri = np.zeros([len(no_vertices), len(shape)])
d_fmm_tri = np.zeros([len(no_vertices), len(shape)])

t_tri_split_gen = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_split = np.zeros([len(no_vertices), len(shape)])
t_fmm_split = np.zeros([len(no_vertices), len(shape)])
d_dijkstra_split = np.zeros([len(no_vertices), len(shape)])
d_fmm_split = np.zeros([len(no_vertices), len(shape)])

# Compute the test example distance for each resolution, shape and mesh type.
for i in range(len(no_vertices)):
    
    # Generate icosahedral triangulation on the unit sphere
    tic = time.perf_counter()
    vertex_ico = leod.triangulation_sphere.triangulate_sphere(1.0, no_divisions[i])
    toc = time.perf_counter()
    t_tri_gen[i] = toc - tic
    
    for j in range(len(shape)):
        
        ### Polar grid with 4-neighbours
        # Generate 4-neighbour polar grid.
        tic = time.perf_counter()
        vertex, polar_grid = leod.fmm_polar_graph.generate_polar_graph(shape[j], no_theta[i], no_phi[i], is_connect_8=False)
        # Precalculation on 4-neighbour polar grid.
        leod.fmm_precalculation.precalculate_grid(vertex)
        toc = time.perf_counter()
        t_polar_gen_4[i][j] = toc - tic
        # Dijkstra's algorithm on the 4-neighbour grid.
        tic = time.perf_counter()
        d_dijkstra_4[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                           order=0, graph_type="polar", grid=polar_grid, is_radians=False)
        toc = time.perf_counter()
        t_dijkstra_4[i][j] = toc - tic
        # Fast marching on the 4-neighbour grid.
        tic = time.perf_counter()
        d_fmm_4[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                      order=1, graph_type="polar", grid=polar_grid, is_radians=False)
        toc = time.perf_counter()
        t_fmm_4[i][j] = toc - tic
        
        
        ### Polar grid with 8-neighbours
        # Generate 8-neighbour polar grid.
        tic = time.perf_counter()
        vertex, polar_grid = leod.fmm_polar_graph.generate_polar_graph(shape[j], no_theta[i], no_phi[i], is_connect_8=True)
        # Precalculation on 8-neighbour polar grid.
        leod.fmm_precalculation.precalculate_grid(vertex)
        toc = time.perf_counter()
        t_polar_gen_8[i][j] = toc - tic
        # Dijkstra's algorithm on the 8-neighbour grid.
        tic = time.perf_counter()
        d_dijkstra_8[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                           order=0, graph_type="polar", grid=polar_grid, is_radians=False)
        toc = time.perf_counter()
        t_dijkstra_8[i][j] = toc - tic
        # Fast marching on the 8-neighbour grid.
        tic = time.perf_counter()
        d_fmm_8[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                      order=1, graph_type="polar", grid=polar_grid, is_radians=False)
        toc = time.perf_counter()
        t_fmm_8[i][j] = toc - tic
        
        ### Icosahedral triangulation
        vertex = copy.deepcopy(vertex_ico)
        # Scale to ellipsoid
        tic = time.perf_counter()
        for k in range(len(vertex)):
            vertex[k].carts[0] *= shape[j].a_axis
            vertex[k].carts[1] *= shape[j].b_axis
            vertex[k].carts[2] *= shape[j].c_axis
        # Precalculation on triangulation
        leod.fmm_precalculation.precalculate_grid(vertex)
        toc = time.perf_counter()
        t_tri_scale_gen[i][j] = t_tri_gen[i] + toc - tic
        # Dijkstra's algorithm on triangulation.
        tic = time.perf_counter()
        d_dijkstra_tri[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                             order=0, graph_type="tri", grid=-1, is_radians=False)
        toc = time.perf_counter()
        t_dijkstra_tri[i][j] = toc - tic
        # Fast marching on triangulation.
        tic = time.perf_counter()
        d_fmm_tri[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                        order=1, graph_type="tri", grid=-1, is_radians=False)
        toc = time.perf_counter()
        t_fmm_tri[i][j] = toc - tic
        
        ### Icosahedral triangulation with splitting
        # Split obtuse angles.
        tic = time.perf_counter()
        leod.fmm_precalculation.split_obtuse_angles(vertex)
        toc = time.perf_counter()
        t_tri_split_gen[i][j] = t_tri_scale_gen[i][j] + toc - tic
        # Dijkstra's algorithm on triangulation.
        tic = time.perf_counter()
        d_dijkstra_split[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                               order=0, graph_type="tri", grid=-1, is_radians=False)
        toc = time.perf_counter()
        t_dijkstra_split[i][j] = toc - tic
        # Fast marching on triangulation.
        tic = time.perf_counter()
        d_fmm_split[i][j], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], vertex, start_point, end_point,
                                                                          order=1, graph_type="tri", grid=-1, is_radians=False)
        toc = time.perf_counter()
        t_fmm_split[i][j] = toc - tic
        
        
# Write distances and run-times to files.
# Use one file per shape.
for j in range(len(shape)):
    file_name = 'data/resolution/resolution_' + shape_name[j] + '.txt'
    with open(file_name, mode="w", newline='') as f:
        
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        top_row = ['d_Dijkstra_4', 'd_FMM_4', 'd_Dijkstra_8', 'd_FMM_8',
                   'd_Dijkstra_Tri', 'd_FMM_Tri', 'd_Dijkstra_Tri_split', 'd_FMM_Tri_split',
                   't_Gen_4', 't_Dijkstra_4', 't_FMM_4', 't_Gen_8', 't_Dijkstra_8', 't_FMM_8',
                   't_Gen_Tri', 't_Dijkstra_Tri', 't_FMM_Tri',
                   't_Gen_Tri_Split', 't_Dijkstra_Tri_Split', 't_FMM_Tri_Split',]
        writer.writerow(top_row)
        
        for i in range(len(no_vertices)):
            row = [ d_dijkstra_4[i][j], d_fmm_4[i][j], d_dijkstra_8[i][j], d_fmm_8[i][j],
                    d_dijkstra_tri[i][j], d_fmm_tri[i][j], d_dijkstra_split[i][j], d_fmm_split[i][j],
                    t_polar_gen_4[i][j], t_dijkstra_4[i][j], t_fmm_4[i][j],
                    t_polar_gen_8[i][j], t_dijkstra_8[i][j], t_fmm_8[i][j],
                    t_tri_scale_gen[i][j], t_dijkstra_tri[i][j], t_fmm_tri[i][j],
                    t_tri_split_gen[i][j], t_dijkstra_split[i][j], t_fmm_split[i][j]]
            writer.writerow(row)        

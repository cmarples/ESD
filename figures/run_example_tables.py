# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 18:27:10 2022

@author: Callum Marples

Calculate distances for a set of examples.

Start and end point angles are given in degrees

Five different shapes are used and for each, an appropriate alternate method
is used for comparison.

This generates data for Tables 1-5 of the Geodesic Paper.
"""

import leod
import math
import numpy as np
import copy
import csv
import sys
import os

save_flag = True

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

### Mesh resolution
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

### Define start and end points
start_point = [ [50.0, 60.0],
                [90.0, 0.0],
                [90.0, 0.0],
                [10.0, 25.0],
                [10.0, 25.0],
                [120.0, 40.0],
                [70.0, 80.0],
                [0.0, 0.0],
                [0.0, 0.0],
                [0.0, 0.0] ]

end_point = [ [90.0, 0.0],
                [50.0, 60.0],
                [130.0, 60.0],
                [20.0, 35.0],
                [10.0, 190.0],
                [60.0, 220.0],
                [40.0, 340.0],
                [180.0, 0.0],
                [90.0, 0.0],
                [90.0, 90.0] ]

### Initialise lists
distance_fmm = np.zeros([5, 4, 10])
# 5 shapes, 4 meshes (polar4, polar8, ico, ico split) and 10 examples.
d_sphere = np.zeros([2, 10])   # True and GeographicLib.
d_wgs84 = np.zeros([2, 10])    # Boundary value method and GeographicLib.
d_oblate = np.zeros([2, 10])   # Boundary value method and GeographicLib.
d_prolate = np.zeros([2, 10])  # Boundary value method and GeographicLib.
d_triaxial = np.zeros([10])    # Boundary value method only.

### Generate and precalculate meshes
polar_4_vertices = [0] * 5
polar_8_vertices = [0] * 5
polar_4_grid = [0] * 5
polar_8_grid = [0] * 5
ico_vertices = [0] * 5
ico_split_vertices = [0] * 5
vertex_sphere_ico = leod.triangulation_sphere.triangulate_sphere(1.0, no_divisions)
for j in range(len(shape)):
    # Generate
    polar_4_vertices[j], polar_4_grid[j] = leod.fmm_polar_graph.generate_polar_graph(shape[j], no_theta, no_phi, is_connect_8=False)
    polar_8_vertices[j], polar_8_grid[j] = leod.fmm_polar_graph.generate_polar_graph(shape[j], no_theta, no_phi, is_connect_8=True)
    ico_vertices[j] = copy.deepcopy(vertex_sphere_ico)
    for k in range(len(ico_vertices[j])):
            ico_vertices[j][k].carts[0] *= shape[j].a_axis
            ico_vertices[j][k].carts[1] *= shape[j].b_axis
            ico_vertices[j][k].carts[2] *= shape[j].c_axis
    ico_split_vertices[j] = copy.deepcopy(ico_vertices[j])
    # Precalculate
    leod.fmm_precalculation.precalculate_grid(polar_4_vertices[j])
    leod.fmm_precalculation.precalculate_grid(polar_8_vertices[j])
    leod.fmm_precalculation.precalculate_grid(ico_vertices[j])
    leod.fmm_precalculation.precalculate_grid(ico_split_vertices[j])
    leod.fmm_precalculation.split_obtuse_angles(ico_split_vertices[j])
print('Meshes generated')

### Calculate fast marching distances
for j in range(len(shape)):
    for i in range(len(start_point)):
        distance_fmm[j][0][i], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], polar_4_vertices[j], start_point[i], end_point[i],
                                                                              order=1, graph_type="polar", grid=polar_4_grid[j], is_radians=False)
        distance_fmm[j][1][i], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], polar_8_vertices[j], start_point[i], end_point[i],
                                                                              order=1, graph_type="polar", grid=polar_8_grid[j], is_radians=False)
        distance_fmm[j][2][i], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], ico_vertices[j], start_point[i], end_point[i],
                                                                              order=1, graph_type="tri", grid=-1, is_radians=False)
        distance_fmm[j][3][i], fmm = leod.fmm_callers.calculate_pair_distance(shape[j], ico_split_vertices[j], start_point[i], end_point[i],
                                                                              order=1, graph_type="tri", grid=-1, is_radians=False)
print('FMM distances calculated')     
 
### Calculate alternate method distances
conv = math.pi / 180.0 # Degrees to radians.
for i in range(len(start_point)):
    # True (sphere only)
    d_sphere[0][i] = leod.sphere_geodesics.great_circle_distance(1.0, start_point[i][0]*conv, start_point[i][1]*conv,
                                                                 end_point[i][0]*conv, end_point[i][1]*conv)
    # GeographicLib (sphere and spheroids).
    if geo_flag == True:
        d_sphere[1][i] = leod.spheroid_geodesics.spheroid_geo_distance(shape[0], start_point[i][0]*conv, start_point[i][1]*conv,
                                                                       end_point[i][0]*conv, end_point[i][1]*conv)
        d_wgs84[1][i] = leod.spheroid_geodesics.spheroid_geo_distance(shape[1], start_point[i][0]*conv, start_point[i][1]*conv,
                                                                       end_point[i][0]*conv, end_point[i][1]*conv)
        d_oblate[1][i] = leod.spheroid_geodesics.spheroid_geo_distance(shape[2], start_point[i][0]*conv, start_point[i][1]*conv,
                                                                       end_point[i][0]*conv, end_point[i][1]*conv)
        d_prolate[1][i] = leod.spheroid_geodesics.spheroid_geo_distance(shape[3], start_point[i][0]*conv, start_point[i][1]*conv,
                                                                       end_point[i][0]*conv, end_point[i][1]*conv)
        
    # Boundary value method (spheroid and triaxial).
    try:
        d = leod.triaxial_geodesics.boundary_value_method(shape[1], start_point[i][0]*conv, start_point[i][1]*conv,
                                                          end_point[i][0]*conv, end_point[i][1]*conv,
                                                          tol=1e-12, Jacobi=False, n = 20000)
        d_wgs84[0][i] = d[0]
    except:
        d_wgs84[0][i] = -1
        
    try:
        d = leod.triaxial_geodesics.boundary_value_method(shape[2], start_point[i][0]*conv, start_point[i][1]*conv,
                                                          end_point[i][0]*conv, end_point[i][1]*conv,
                                                          tol=1e-12, Jacobi=False, n = 20000)
        d_oblate[0][i] = d[0]
    except:
        d_oblate[0][i] = -1
        
    try:
        d = leod.triaxial_geodesics.boundary_value_method(shape[3], start_point[i][0]*conv, start_point[i][1]*conv,
                                                          end_point[i][0]*conv, end_point[i][1]*conv,
                                                          tol=1e-12, Jacobi=False, n = 20000)
        d_prolate[0][i] = d[0]
    except:
        d_prolate[0][i] = -1
        
    try:
        d = leod.triaxial_geodesics.boundary_value_method(shape[4], start_point[i][0]*conv, start_point[i][1]*conv,
                                                          end_point[i][0]*conv, end_point[i][1]*conv,
                                                          tol=1e-12, Jacobi=False, n = 20000)
        d_triaxial[i] = d[0]
    except:
        d_triaxial[i] = -1
print('Alternate distances generated')        
        
### Write to files
if save_flag == True:
    # Sphere
    file_name1 = 'data/resolution/table_sphere.csv'
    with open(file_name1, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['$theta_0$', '$phi_0$', '$theta_1$', '$phi_1$', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'True', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [start_point[i][0], start_point[i][1], end_point[i][0], end_point[i][1],
                   distance_fmm[0][0][i], distance_fmm[0][1][i], distance_fmm[0][2][i], distance_fmm[0][3][i],
                   d_sphere[0][i], d_sphere[1][i]]
            writer.writerow(row)
    # WGS84
    file_name2 = 'data/resolution/table_wgs84.csv'
    with open(file_name2, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['$theta_0$', '$phi_0$', '$theta_1$', '$phi_1$', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [start_point[i][0], start_point[i][1], end_point[i][0], end_point[i][1],
                   distance_fmm[1][0][i], distance_fmm[1][1][i], distance_fmm[1][2][i], distance_fmm[1][3][i],
                   d_wgs84[0][i], d_wgs84[1][i]]
            writer.writerow(row)
    # Oblate
    file_name3 = 'data/resolution/table_oblate.csv'
    with open(file_name3, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['$theta_0$', '$phi_0$', '$theta_1$', '$phi_1$', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [start_point[i][0], start_point[i][1], end_point[i][0], end_point[i][1],
                   distance_fmm[2][0][i], distance_fmm[2][1][i], distance_fmm[2][2][i], distance_fmm[2][3][i],
                   d_oblate[0][i], d_oblate[1][i]]
            writer.writerow(row)
    # Prolate
    file_name3 = 'data/resolution/table_prolate.csv'
    with open(file_name3, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['$theta_0$', '$phi_0$', '$theta_1$', '$phi_1$', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [start_point[i][0], start_point[i][1], end_point[i][0], end_point[i][1],
                   distance_fmm[3][0][i], distance_fmm[3][1][i], distance_fmm[3][2][i], distance_fmm[3][3][i],
                   d_prolate[0][i], d_prolate[1][i]]
            writer.writerow(row)
    # Triaxial
    file_name3 = 'data/resolution/table_triaxial.csv'
    with open(file_name3, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['$theta_0$', '$phi_0$', '$theta_1$', '$phi_1$', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [start_point[i][0], start_point[i][1], end_point[i][0], end_point[i][1],
                   distance_fmm[4][0][i], distance_fmm[4][1][i], distance_fmm[4][2][i], distance_fmm[4][3][i],
                   d_triaxial[i]]
            writer.writerow(row)

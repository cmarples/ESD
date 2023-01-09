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

modulename = 'leod.geo.spheroid'
if modulename in sys.modules:
    geo_flag = True
else:
    geo_flag = False

os.chdir("..")

### Define shapes.
shape = [ leod.shape.EllipsoidShape(1.0, 1.0, 1.0),   # Sphere.
          leod.shape.EllipsoidShape(6378137.0, 6378137.0, 6356752.3142),  # Earth (WGS84 ellipsoid).
          leod.shape.EllipsoidShape(2.0, 2.0, 1.0),   # Oblate.
          leod.shape.EllipsoidShape(1.0, 1.0, 2.0),   # Prolate.
          leod.shape.EllipsoidShape(3.0, 2.0, 1.0) ]  # Triaxial.
shape_name = ['sphere', 'wgs84', 'oblate', 'prolate', 'triaxial']

# Normalise shapes so that (abc) = 1.
# Unit sphere already normalised, do not want to do this for the WGS84 example.
shape[2].normalise()
shape[3].normalise()
shape[4].normalise()

conv = math.pi / 180.0 # Degrees to radians.

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
start_point = np.array([ [50.0, 60.0],
                         [90.0, 0.0],
                         [90.0, 0.0],
                         [10.0, 25.0],
                         [10.0, 25.0],
                         [120.0, 40.0],
                         [70.0, 80.0],
                         [0.0, 0.0],
                         [0.0, 0.0],
                         [0.0, 0.0],
                         [38.5, 359,9],   # London
                         [54.3, 139.8] ]) # Tokyo

end_point = np.array([ [90.0, 0.0],
                       [50.0, 60.0],
                       [130.0, 60.0],
                       [20.0, 35.0],
                       [10.0, 190.0],
                       [60.0, 220.0],
                       [40.0, 340.0],
                       [180.0, 0.0],
                       [90.0, 0.0],
                       [90.0, 90.0],
                       [112.9, 316.8],  # Rio de Janeiro
                       [41.1, 2.3] ])   # Paris

n = len(start_point)
for i in range(n):
    for j in range(2):
        start_point[i][j] *= conv
        end_point[i][j] *= conv

### Initialise lists
distance_fmm = np.zeros([5, 4, n])
# 5 shapes, 4 meshes (polar4, polar8, ico, ico split) and 12 examples.
d_sphere = np.zeros([2, n])   # True and GeographicLib.
d_wgs84 = np.zeros([2, n])    # Boundary value method and GeographicLib.
d_oblate = np.zeros([2, n])   # Boundary value method and GeographicLib.
d_prolate = np.zeros([2, n])  # Boundary value method and GeographicLib.
d_triaxial = np.zeros([n])    # Boundary value method only.

### Generate and precalculate meshes
polar_4_grid = [0] * 5
polar_8_grid = [0] * 5
ico_grid = [0] * 5
ico_split_grid = [0] * 5

for j in range(len(shape)):
    # Generate
    polar_4_grid[j] = leod.fmm.grid_pol.gen_pol_grid(no_theta, no_phi, shape[j], is_connect_8=False)
    polar_8_grid[j] = leod.fmm.grid_pol.gen_pol_grid(no_theta, no_phi, shape[j], is_connect_8=True)
    ico_grid[j] = leod.fmm.grid_ico.gen_ico_grid(no_divisions, shape[j], is_generic=False)
    ico_split_grid[j] = leod.fmm.grid_ico.gen_ico_grid(no_divisions, shape[j], is_generic=False)
    leod.fmm.grid_pre.split_obtuse_angles(ico_split_grid[j].vertex)
print('Grids generated')


### Calculate fast marching distances
for j in range(len(shape)):
    for i in range(len(start_point)):
        distance_fmm[j][0][i], fmm = leod.fmm.callers.distance_pair(shape[j], polar_4_grid[j], start_point[i], end_point[i], is_dijkstra=False, is_radians=True)
        distance_fmm[j][1][i], fmm = leod.fmm.callers.distance_pair(shape[j], polar_8_grid[j], start_point[i], end_point[i], is_dijkstra=False, is_radians=True)
        distance_fmm[j][2][i], fmm = leod.fmm.callers.distance_pair(shape[j], ico_grid[j], start_point[i], end_point[i], is_dijkstra=False, is_radians=True)
        distance_fmm[j][3][i], fmm = leod.fmm.callers.distance_pair(shape[j], ico_split_grid[j], start_point[i], end_point[i], is_dijkstra=False, is_radians=True)
print('FMM distances calculated')     
 
### Calculate alternate method distances

for i in range(len(start_point)):
    # True (sphere only)
    d_sphere[0][i] = leod.geo.sphere.gc_dist(shape[0].a_axis, start_point[i], end_point[i])

    # GeographicLib (sphere and spheroids).
    if geo_flag == True:
        
        d_sphere[1][i] = leod.geo.spheroid.geo_dist(shape[0], start_point[i], end_point[i])
        d_wgs84[1][i] = leod.geo.spheroid.geo_dist(shape[1], start_point[i], end_point[i])
        d_oblate[1][i] = leod.geo.spheroid.geo_dist(shape[2], start_point[i], end_point[i])
        d_prolate[1][i] = leod.geo.spheroid.geo_dist(shape[3], start_point[i], end_point[i])
        
    # Boundary value method (spheroid and triaxial).
    try:
        d = leod.geo.triaxial.bvm_dist(shape[1], start_point[i], end_point[i], tol=1e-12, Jacobi=False, n = 20000)
        d_wgs84[0][i] = d[0]
    except:
        d_wgs84[0][i] = -1
        
    try:
        d = leod.geo.triaxial.bvm_dist(shape[2], start_point[i], end_point[i], tol=1e-12, Jacobi=False, n = 20000)
        d_oblate[0][i] = d[0]
    except:
        d_oblate[0][i] = -1
        
    try:
        d = leod.geo.triaxial.bvm_dist(shape[3], start_point[i], end_point[i], tol=1e-12, Jacobi=False, n = 20000)
        d_prolate[0][i] = d[0]
    except:
        d_prolate[0][i] = -1
        
    try:
        d = leod.geo.triaxial.bvm_dist(shape[4], start_point[i], end_point[i], tol=1e-12, Jacobi=False, n = 20000)
        d_triaxial[i] = d[0]
    except:
        d_triaxial[i] = -1
print('Alternate distances generated')        
        
### Write to files
# Get start and end points back in degrees
start_point = np.array([ [50.0, 60.0],
                         [90.0, 0.0],
                         [90.0, 0.0],
                         [10.0, 25.0],
                         [10.0, 25.0],
                         [120.0, 40.0],
                         [70.0, 80.0],
                         [0.0, 0.0],
                         [0.0, 0.0],
                         [0.0, 0.0],
                         [38.5, 359,9],   # London
                         [54.3, 139.8] ]) # Tokyo

end_point = np.array([ [90.0, 0.0],
                       [50.0, 60.0],
                       [130.0, 60.0],
                       [20.0, 35.0],
                       [10.0, 190.0],
                       [60.0, 220.0],
                       [40.0, 340.0],
                       [180.0, 0.0],
                       [90.0, 0.0],
                       [90.0, 90.0],
                       [112.9, 316.8],  # Rio de Janeiro
                       [41.1, 2.3] ])   # Paris

if save_flag == True:
    # Sphere
    file_name1 = 'data/geodesics/tables/table_sphere.csv'
    with open(file_name1, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['(theta_0, phi_0)', '(theta_1, phi_1)', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'True', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [(start_point[i][0], start_point[i][1]), (end_point[i][0], end_point[i][1]),
                   distance_fmm[0][0][i], distance_fmm[0][1][i], distance_fmm[0][2][i], distance_fmm[0][3][i],
                   d_sphere[0][i], d_sphere[1][i]]
            writer.writerow(row)
    # WGS84
    file_name2 = 'data/geodesics/tables/table_wgs84.csv'
    with open(file_name2, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['(theta_0, phi_0)', '(theta_1, phi_1)', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [(start_point[i][0], start_point[i][1]), (end_point[i][0], end_point[i][1]),
                   distance_fmm[1][0][i], distance_fmm[1][1][i], distance_fmm[1][2][i], distance_fmm[1][3][i],
                   d_wgs84[0][i], d_wgs84[1][i]]
            writer.writerow(row)
    # Oblate
    file_name3 = 'data/geodesics/tables/table_oblate.csv'
    with open(file_name3, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['(theta_0, phi_0)', '(theta_1, phi_1)', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [(start_point[i][0], start_point[i][1]), (end_point[i][0], end_point[i][1]),
                   distance_fmm[2][0][i], distance_fmm[2][1][i], distance_fmm[2][2][i], distance_fmm[2][3][i],
                   d_oblate[0][i], d_oblate[1][i]]
            writer.writerow(row)
    # Prolate
    file_name4 = 'data/geodesics/tables/table_prolate.csv'
    with open(file_name4, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['(theta_0, phi_0)', '(theta_1, phi_1)', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [(start_point[i][0], start_point[i][1]), (end_point[i][0], end_point[i][1]),
                   distance_fmm[3][0][i], distance_fmm[3][1][i], distance_fmm[3][2][i], distance_fmm[3][3][i],
                   d_prolate[0][i], d_prolate[1][i]]
            writer.writerow(row)
    # Triaxial
    file_name5 = 'data/geodesics/tables/table_triaxial.csv'
    with open(file_name5, mode="w", newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        top_row = ['(theta_0, phi_0)', '(theta_1, phi_1)', 'Polar 4', 'Polar 8', 'Icosahedral', 'Icosahedral (Split)', 'BVM', 'GeographicLib']
        writer.writerow(top_row)
        for i in range(len(start_point)):
            row = [(start_point[i][0], start_point[i][1]), (end_point[i][0], end_point[i][1]),
                   distance_fmm[4][0][i], distance_fmm[4][1][i], distance_fmm[4][2][i], distance_fmm[4][3][i],
                   d_triaxial[i]]
            writer.writerow(row)
            
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 14:46:25 2022

@author: Callum Marples

Calculate distance for a particular example, at a range of resolutions.

Example: (theta, phi) = (50, 60) -> (90, 0)
(above angles in degrees)

Five different shapes are used and for each, an appropriate alternate method
is used for comparison.
"""

import leod
import math
import time
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

# Define shapes.
shape = [ leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 1.0),   # Sphere.
          leod.ellipsoid_shape.EllipsoidShape(6378137.0, 6378137.0, 6356752.3142),  # Earth (WGS84 ellipsoid).
          leod.ellipsoid_shape.EllipsoidShape(2.0, 2.0, 1.0),   # Oblate.
          leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 2.0),   # Prolate.
          leod.ellipsoid_shape.EllipsoidShape(3.0, 2.0, 1.0) ]  # Triaxial.
shape_name = ['Sphere', 'WGS84', 'Oblate', 'Prolate', 'Triaxial']

# Normalise shapes so that (abc) = 1.
# Unit sphere already normalised, do not want to do this for the WGS84 example.
shape[2].normalise()
shape[3].normalise()
shape[4].normalise()

# Define start and end points of the geodesic example.
start_point = [50.0, 60.0]
end_point = [90.0, 0.0]
conv = math.pi / 180.0 # Degrees to radians.

dist = []
taxicab = []
run_time = []
# Compute the distance on each shape in turn, for each alternate method.
for i in range(len(shape)):
    dist.append([])
    run_time.append([])
    d_true = 'N/A'
    d_geo = 'N/A'
    d_bvm = 'N/A'
    d_taxi = 'N/A'
    t_true = 'N/A'
    t_geo = 'N/A'
    t_bvm = 'N/A'
    t_taxi = 'N/A'
    if shape[i].is_sphere() == True:
        # Sphere (can calculate true value and call GeographicLib).
        tic = time.perf_counter()
        d_true = leod.sphere_geodesics.great_circle_distance(shape[i].a_axis, start_point[0]*conv, start_point[1]*conv,
                                       end_point[0]*conv, end_point[1]*conv)
        toc = time.perf_counter()
        t_true = toc - tic
        # Sphere taxicab distance
        tic = time.perf_counter()
        d_taxi = leod.taxicab_distance.taxicab_sphere(1.0, start_point[0]*conv, start_point[1]*conv,
                                       end_point[0]*conv, end_point[1]*conv)
        toc = time.perf_counter()
        t_taxi = toc - tic
        if geo_flag == True:
            tic = time.perf_counter()
            d_geo = leod.spheroid_geodesics.spheroid_geo_distance(shape[i], start_point[0]*conv, start_point[1]*conv,
                                                                  end_point[0]*conv, end_point[1]*conv)
            toc = time.perf_counter()
            t_geo = toc - tic
    else:    
        if shape[i].is_spheroid() == True:
            # Spheroid (GeographicLib and Boundary Value Method).
            if geo_flag == True:
                tic = time.perf_counter()
                d_geo = leod.spheroid_geodesics.spheroid_geo_distance(shape[i], start_point[0]*conv, start_point[1]*conv,
                                                                      end_point[0]*conv, end_point[1]*conv)
                toc = time.perf_counter()
                t_geo = toc - tic
            # Spheroid taxicab distance
            tic = time.perf_counter()
            d_taxi = leod.taxicab_distance.taxicab_spheroid(shape[i].a_axis, shape[i].c_axis, start_point[0]*conv, start_point[1]*conv,
                                                            end_point[0]*conv, end_point[1]*conv)
            toc = time.perf_counter()
            t_taxi = toc - tic
        # Triaxial (Boundary Value Method only).
        tic = time.perf_counter()
        d_bvm = leod.triaxial_geodesics.boundary_value_method(shape[i], start_point[0]*conv, start_point[1]*conv,
                                                              end_point[0]*conv, end_point[1]*conv,
                                                              tol=1e-12, Jacobi=False, n = 20000)
        d_bvm = d_bvm[0]
        toc = time.perf_counter()
        t_bvm = toc - tic
    dist[i] = [d_true, d_geo, d_bvm, d_taxi]
    run_time[i] = [t_true, t_geo, t_bvm, t_taxi]
    
                
# Write alternate method results and shape information to a file
if save_flag == True:
    file_name = 'data/resolution/resolution_comparison.txt'
    with open(file_name, mode="w", newline='') as f:
        
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        top_row = ['Shape', 'a', 'b', 'c', 'True', 'GeographicLib', 'Boundary Value', 't_True', 't_Geo', 't_BV']
        writer.writerow(top_row)
        
        for i in range(len(shape)):
            row = [shape_name[i], shape[i].a_axis, shape[i].b_axis, shape[i].c_axis, dist[i][0], dist[i][1], dist[i][2], dist[i][3], run_time[i][0], run_time[i][1], run_time[i][2], run_time[i][3]]
            writer.writerow(row)
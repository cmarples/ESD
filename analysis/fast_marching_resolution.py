# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 15:20:10 2022

@author: Callum Marples
"""

import math
import numpy as np
import csv
import time

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid
from leod.geo_fmm import GeoFMM
from leod.sphere_geodesics import great_circle_distance
from leod.spheroid_geodesics import spheroid_geo_distance
from leod.triaxial_geodesics import boundary_value_method

# Get shape type
print('Distance as a function of resolution.')
print('Start point: (theta, phi) = (90, 0)')
print('End point  : (theta, phi) = (50, 60)')
print('The following shapes are used:')
print('')
print('    Sphere       (1.0, 1.0, 1.0)')
print('    Prolate      (1.0, 1.0, 2.0) / norm')
print('    Triaxial     (3.0, 2.0, 1.0) / norm')
print('    Earth WGS84  (6378137, 6378137, 6356752.3142)')

# Start and end points
th_0 = 90.0 * math.pi / 180.0
ph_0 = 0.0
th_1 = 50.0 * math.pi / 180.0
ph_1 = 60.0 * math.pi / 180.0
    
for shape_count in range(4):
        
    # File name, ellipsoid axes and EllipsoidShape object        
    if shape_count == 0:
        E = EllipsoidShape(1.0, 1.0, 1.0)
        file_name = "data/distance_vs_resolution_sphere.csv"
        s = great_circle_distance(1.0, th_0, ph_0, th_1, ph_1)
    elif shape_count == 1:
        E = EllipsoidShape(1.0, 1.0, 2.0)
        file_name = "data/distance_vs_resolution_prolate.csv"
        E.normalise()
        s = spheroid_geo_distance(E, th_0, ph_0, th_1, ph_1)
    elif shape_count == 2:
        E = EllipsoidShape(3.0, 2.0, 1.0)
        file_name = "data/distance_vs_resolution_triaxial.csv"
        E.normalise()
        s = boundary_value_method(E, th_0, ph_0, th_1, ph_1)
    else:
        E = EllipsoidShape(6378137, 6378137, 6356752.3142)
        file_name = "data/distance_vs_resolution_earth.csv"
        s = spheroid_geo_distance(E, th_0, ph_0, th_1, ph_1)
        
    n_res = np.arange(25, 525, step=25)
    m = len(n_res)
    d_1_list = [0] * m
    t_1_list = [0] * m
    d_2_list = [0] * m
    t_2_list = [0] * m
    d_3_list = [0] * m
    d_4_list = [0] * m
    t_3_list = [0] * m
    t_4_list = [0] * m
    
    x_ref = 10
    n_div = 5
    
    for i in range(m):
        
        # First order FMM #
        G1 = GeoGrid(E, n_res[i], n_res[i]+1)
        F1 = GeoFMM(G1, th_0, ph_0)
        # Compute geodesic
        tic = time.perf_counter()
        d_1_list[i] = F1.calculate_geodesics(1, th_1, ph_1)
        toc = time.perf_counter()
        t_1_list[i] = toc - tic
        
        # Second order FMM #
        G2 = GeoGrid(E, n_res[i], n_res[i]+1)
        F2 = GeoFMM(G2, th_0, ph_0)
        # Compute geodesic
        tic = time.perf_counter()
        d_2_list[i] = F2.calculate_geodesics(2, th_1, ph_1)
        toc = time.perf_counter()
        t_2_list[i] = toc - tic
        
        # First order refined FMM #
        G3 = GeoGrid(E, n_res[i], n_res[i]+1)
        F3 = GeoFMM(G3, th_0, ph_0)
        # Compute geodesic
        tic = time.perf_counter()
        d_3_list[i] = F3.calculate_geodesics(1, th_1, ph_1, True, x_ref, n_div, n_div)
        toc = time.perf_counter()
        t_3_list[i] = toc - tic
        
        # Second order refined FMM #
        G4 = GeoGrid(E, n_res[i], n_res[i]+1)
        F4 = GeoFMM(G4, th_0, ph_0)
        # Compute geodesic
        tic = time.perf_counter()
        d_4_list[i] = F4.calculate_geodesics(2, th_1, ph_1, True, x_ref, n_div, n_div)
        toc = time.perf_counter()
        t_4_list[i] = toc - tic
           
    with open(file_name, mode="w", newline='') as f:
        
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        writer.writerow(d_1_list)
        writer.writerow(d_2_list)
        writer.writerow(d_3_list)
        writer.writerow(d_4_list)
        
        writer.writerow(t_1_list)
        writer.writerow(t_2_list)
        writer.writerow(t_3_list)
        writer.writerow(t_4_list)
            
        writer.writerow([E.a_axis, E.b_axis, E.c_axis])
        writer.writerow([s])
        
    print("")    
    if shape_count == 0:
        print("Sphere completed")
    elif shape_count == 1:
        print("Prolate completed")
    elif shape_count == 2:
        print("Triaxial completed")
    else:
        print("Earth completed")
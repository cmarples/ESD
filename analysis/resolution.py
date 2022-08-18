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

# Get shape type
print('')
print('Enter one of hte following for shape type:')
print('')
print('    "s" : sphere   (1.0, 1.0, 1.0)')
print('    "p" : prolate  (1.0, 1.0, 2.0) / norm')
print('    "t" : triaxial (3.0, 2.0, 1.0) / norm')
print('    "e" : earth    (6378137, 6378137, 6356752.3142)')
shape_type = input()
if shape_type != "s" and shape_type != "p" and shape_type != "t" and shape_type != "e":
    while shape_type != "s" and shape_type != "p" and shape_type != "t" and shape_type != "e":
        print('')
        print('Invalid shape type. Select "s", "p", "t", or "e":')
        print('')
        shape_type = input()
        
# File name, ellipsoid axes and EllipsoidShape object        
file_name = "data/" + shape_type + "_fmm.txt"
if shape_type == "s":
    E = EllipsoidShape(1.0, 1.0, 1.0)
elif shape_type == "p":
    E = EllipsoidShape(1.0, 2.0, 2.0, True)
elif shape_type == "t":
    E = EllipsoidShape(3.0, 2.0, 1.0, True)
else:
    E = EllipsoidShape(6378137, 6378137, 6356752.3142)

# Start and end points
th_0 = 90.0 * math.pi / 180.0
ph_0 = 0.0
th_1 = 50.0 * math.pi / 180.0
ph_1 = 60.0 * math.pi / 180.0

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
    G1 = GeoGrid(E, n_res[i], n_res[i])
    F1 = GeoFMM(G1, th_0, ph_0)
    # Compute geodesic
    tic = time.perf_counter()
    d_1_list[i] = F1.calculate_geodesics(1, th_1, ph_1)
    toc = time.perf_counter()
    t_1_list = toc - tic
    
    # Second order FMM #
    G2 = GeoGrid(E, n_res[i], n_res[i])
    F2 = GeoFMM(G2, th_0, ph_0)
    # Compute geodesic
    tic = time.perf_counter()
    d_2_list[i] = F2.calculate_geodesics(2, th_1, ph_1)
    toc = time.perf_counter()
    t_2_list = toc - tic
    
    # First order refined FMM #
    G3 = GeoGrid(E, n_res[i], n_res[i])
    F3 = GeoFMM(G3, th_0, ph_0)
    # Compute geodesic
    tic = time.perf_counter()
    d_3_list[i] = F3.calculate_geodesics(1, th_1, ph_1, True, x_ref, n_div, n_div)
    toc = time.perf_counter()
    t_3_list = toc - tic
    
    # Second order refined FMM #
    G4 = GeoGrid(E, n_res[i], n_res[i])
    F4 = GeoFMM(G4, th_0, ph_0)
    # Compute geodesic
    tic = time.perf_counter()
    d_4_list[i] = F4.calculate_geodesics(2, th_1, ph_1, True, x_ref, n_div, n_div)
    toc = time.perf_counter()
    t_4_list = toc - tic
        
with open(file_name, mode="w") as f:
    
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
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 16:39:45 2022

@author: Callum Marples

This file is used for devopment purposes (to test/debug the code).
It works by generating the necessary objects and calculating the distances
for a simple example.
"""

import os
os.chdir("..")

import math
import time

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid
from leod.geo_fmm import GeoFMM
from leod.sphere_geodesics import great_circle_distance
from leod.spheroid_geodesics import spheroid_geo_distance
from leod.triaxial_geodesics import boundary_value_method
from leod.taxicab_distance import taxicab_distance_sphere

import numpy as np

test_type = 12

if test_type == 1:   # Generate the GeoGrid, GeoPixel and GeoFMM objects
    
    E = EllipsoidShape(3.0, 2.0, 1.0)
    G = GeoGrid(E, 181, 360)
    th = 75.0 * math.pi / 180.0
    ph = 100.0 * math.pi / 180.0
    F = GeoFMM(G, th, ph)
    
elif test_type == 2: # 4-neighbour Dijkstra example
    
    E = EllipsoidShape(3.0, 2.0, 1.0)
    G = GeoGrid(E, 181, 360)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    F = GeoFMM(G, th_0, ph_0)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    # Expect d = 2.8676959922
    d = F.calculate_geodesics(0, th_1, ph_1)
    
elif test_type == 3: # 1st Order FMM example
    
    E = EllipsoidShape(1.0, 1.0, 1.0)
    r = 125
    G = GeoGrid(E, r, r)
    th_0 = 90.0001 * math.pi / 180.0
    ph_0 = 0.0
    F = GeoFMM(G, th_0, ph_0)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    d = F.calculate_geodesics(1, th_1, ph_1)
    
elif test_type == 4: # 2nd Order FMM example
    
    E = EllipsoidShape(1.0, 1.0, 1.0)
    G = GeoGrid(E, 200, 200)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    tic = time.perf_counter()
    F = GeoFMM(G, th_0, ph_0)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    # Expect d = 2.3550639686
    d = F.calculate_geodesics(2, th_1, ph_1)
    toc = time.perf_counter()
    print(toc - tic)
elif test_type == 5: # 2nd Order FMM example (with refinement)
    
    n = 200 # n_theta = n_phi
    x_ref = 10
    no_div = 5
    E = EllipsoidShape(1.0, 1.0, 1.0)
    G = GeoGrid(E, 190, 350)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    tic = time.perf_counter()
    F = GeoFMM(G, th_0, ph_0)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    d = F.calculate_geodesics(2, th_1, ph_1, is_refine=True, refine_range=x_ref, 
                              refine_theta=no_div, refine_phi=no_div)
    toc = time.perf_counter()
    print(toc - tic)
    
    s = great_circle_distance(1.0, th_0, ph_0, th_1, ph_1)
    
elif test_type == 6: # 2nd Order FMM example on a flat grid
    
    n = 200 # n_theta = n_phi
    x_ref = 10
    no_div = 5
    E = EllipsoidShape(1.0, 1.0, 1.0)
    G = GeoGrid(E, 201, 200, is_flat = True)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    tic = time.perf_counter()
    F = GeoFMM(G, th_0, ph_0, is_flat = True)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    d = F.calculate_geodesics(2, th_1, ph_1, is_refine=True, refine_range=x_ref, 
                              refine_theta=no_div, refine_phi=no_div)
    toc = time.perf_counter()
    print(toc - tic)
    
    s = great_circle_distance(1.0, th_0, ph_0, th_1, ph_1)
    
elif test_type == 7: # Sphere geodesic using GeographicLib
    E = EllipsoidShape(1.0, 1.0, 1.0)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    c = great_circle_distance(1.0, th_0, ph_0, th_1, ph_1)
    s = spheroid_geo_distance(E, th_0, ph_0, th_1, ph_1)
    
elif test_type == 8: # Spheroid geodesic using GeographicLib     
    E = EllipsoidShape(6378137.0, 6378137.0, 6356752.3142) # WGS84 ellipsoid
    th_0 = (90.0 - 41.32) * math.pi / 180.0
    ph_0 = 174.81 * math.pi / 180.0
    th_1 = (90.0 + 40.96) * math.pi / 180.0
    ph_1 = -5.50 * math.pi / 180.0
    s = spheroid_geo_distance(E, th_0, ph_0, th_1, ph_1)
    
elif test_type == 9: # Triaxial geodesics using Panou's boundary value method
    E = EllipsoidShape(6378172.0, 6378103.0, 6356753.0) # Earth parameters
    th_0 = 0.0
    ph_0 = 0.0
    th_1 = 0.0
    ph_1 = 90.0 * math.pi / 180.0
    #s = boundary_value_method(E, th_0, ph_0, th_1, ph_1, 1.0e-12, True, 10000)
    #s = s[0]
    conv = math.pi / 180.0
    s = boundary_value_method(E, 30.0*conv, 0.0, -30.0*conv, 175.0*conv, tol=1e-12, Jacobi=True, n = 16000) 

elif test_type == 10: # Prolate and Oblate
   
   n = 200 # n_phi
   th_0 = 90.0 * math.pi / 180.0
   ph_0 = 0.0
   th_1 = 50.0 * math.pi / 180.0
   ph_1 = 60.0 * math.pi / 180.0
   
   # Prolate
   E1 = EllipsoidShape(1.0, 1.0, 2.0)
   E1.normalise()
   G1 = GeoGrid(E1, n+1, n)
   F1 = GeoFMM(G1, th_0, ph_0)
   d1 = F1.calculate_geodesics(2, th_1, ph_1, is_refine=False)
   s1 = spheroid_geo_distance(E1, th_0, ph_0, th_1, ph_1)   
    
   # Oblate
   E2 = EllipsoidShape(2.0, 2.0, 1.0)
   E2.normalise()
   G2 = GeoGrid(E2, n+1, n)
   F2 = GeoFMM(G2, th_0, ph_0)
   d2 = F2.calculate_geodesics(2, th_1, ph_1, is_refine=False)
   s2 = spheroid_geo_distance(E2, th_0, ph_0, th_1, ph_1) 
   
elif test_type == 11: # 8-neighbour Dijkstra
    a = 1.0
    b = 1.0
    c = 1.0
    E = EllipsoidShape(a, b, c)
    G = GeoGrid(E, 200, 200, neighbour8=False)
    
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    F = GeoFMM(G, th_0, ph_0)
    
    d = F.calculate_geodesics(0, th_1, ph_1)
    
    s = taxicab_distance_sphere(a, th_0, ph_0, th_1, ph_1)
    
elif test_type == 12: # Multiple endpoints
    a = 1.0
    b = 1.0
    c = 1.0
    E = EllipsoidShape(a, b, c)
    G = GeoGrid(E, 200, 200, neighbour8=False)
    
    deg2rad = math.pi / 180.0
    th_0 = 90.0 * deg2rad
    ph_0 = 0.0
    th_1 = [50.0 * deg2rad, 100.0 * deg2rad, 50.1 * deg2rad]
    ph_1 = [60.0 * deg2rad, 20.0 * deg2rad, 59.9 * deg2rad]
    
    F = GeoFMM(G, th_0, ph_0)
    
    d = F.calculate_geodesics(2, th_1, ph_1, is_refine=True, refine_range=10, refine_theta=3, refine_phi=3)
    
    s = [-1.0] * 3
    for i in range(3):
        s[i] = taxicab_distance_sphere(a, th_0, ph_0, th_1[i], ph_1[i])
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

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid
from leod.geo_fmm import GeoFMM

test_type = 5

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
    
    E = EllipsoidShape(3.0, 2.0, 1.0)
    G = GeoGrid(E, 181, 360)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    F = GeoFMM(G, th_0, ph_0)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    # Expect d = 2.3854095543
    d = F.calculate_geodesics(1, th_1, ph_1)
    
elif test_type == 4: # 2nd Order FMM example
    
    E = EllipsoidShape(1.0, 1.0, 1.0)
    G = GeoGrid(E, 200, 200)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    F = GeoFMM(G, th_0, ph_0)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    # Expect d = 2.3550639686
    d = F.calculate_geodesics(2, th_1, ph_1)

elif test_type == 5: # 2nd Order FMM example (with refinement)
    
    E = EllipsoidShape(1.0, 1.0, 1.0)
    G = GeoGrid(E, 200, 200)
    th_0 = 90.0 * math.pi / 180.0
    ph_0 = 0.0
    F = GeoFMM(G, th_0, ph_0)
    th_1 = 50.0 * math.pi / 180.0
    ph_1 = 60.0 * math.pi / 180.0
    
    d = F.calculate_geodesics(2, th_1, ph_1, is_refine=True, refine_range=1, 
                              refine_theta=3, refine_phi=3)
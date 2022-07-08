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

test_type = 1

if test_type == 1:   # Generate the GeoGrid, GeoPixel and GeoFMM objects
    
    E = EllipsoidShape(3.0, 2.0, 1.0)
    G = GeoGrid(E, 181, 360)
    th = 75.0 * math.pi / 180.0
    ph = 100.0 * math.pi / 180.0
    F = GeoFMM(G, th, ph)
    

    
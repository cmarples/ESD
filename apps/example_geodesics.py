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

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid

test_type = 1

if test_type == 1: # Generate the GeoGrid and GeoPixel objects
    # Generate shape
    E = EllipsoidShape(3.0, 2.0, 1.0)
    G = GeoGrid(E, 181, 360)
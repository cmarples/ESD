# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:54:00 2022

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
from leod.fmm_vertex import FmmVertex
from leod.fmm_grid import FmmGrid
from leod.fmm_create_polar import generate_polar_grid

test_no = 1

if test_no == 1: # Create ellipsoid grid using polar coordinates
    shape = EllipsoidShape(3.0, 2.0, 1.0)
    grid = generate_polar_grid(shape, 6, 4)
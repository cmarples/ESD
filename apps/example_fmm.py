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
from leod.fmm_polar_graph import generate_polar_graph
from leod.fmm_polar_graph import check_triangles

test_no = 1

if test_no == 1: # Create ellipsoid grid using polar coordinates
    shape = EllipsoidShape(1.0, 1.0, 2.0)
    vertex = generate_polar_graph(shape, 181, 360, is_connect_8=True, is_Dijkstra=False)
    
    # Check that these triangles are acute (and thus valid for the FMM)
    no_obtuse = check_triangles(vertex)
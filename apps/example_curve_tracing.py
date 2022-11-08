# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 15:00:48 2022

@author: Callum Marples

Given ellipsoid (a,b,c) and a sphere of radius r (centred at an ellipsoid
surface point), find a sample of points on the intersecting curve
"""

import os
os.chdir("..")

import math
import time
import numpy as np

from leod.ellipsoid_shape import EllipsoidShape
import leod.tri_curve_tracing as trace

tic = time.perf_counter()
    
r = 0.01
shape = EllipsoidShape(3.0, 2.0, 1.0)

d2r = math.pi / 180.0 # Degrees to radians

th_0 = 90.0 * d2r
ph_0 = 180.0 * d2r
centre = np.array(shape.polar2cart(th_0, ph_0))

th = trace.initial_point(centre, r, shape.a_axis*math.cos(ph_0), shape.b_axis*math.sin(ph_0), 
                         shape.c_axis, math.cos(th_0))

point = np.array(shape.polar2cart(th, ph_0))
x = centre - point
d = math.sqrt(np.dot(x, x))

m = trace.h_implicit(th, ph_0, centre, shape.a_axis, shape.b_axis, shape.c_axis, r*r)

toc = time.perf_counter()
print(toc - tic)



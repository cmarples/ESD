# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 14:32:03 2022

@author: Callum Marples

This file contains code to call routines from the GeographicLib library to
calculate geodesic distnaces on the spheroid.

To run any function in this file, GeographicLib must be installed!
"""

import math
#from .ellipsoid_shape import EllipsoidShape
from geographiclib.geodesic import Geodesic


# Calculate shortest distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a spheroid of axis lengths (a, a, c), as given in an
# EllipsoidShape object
def spheroid_geo_distance(shape, theta_0, phi_0, theta_1, phi_1):
    flattening = (shape.a_axis - shape.c_axis) / shape.a_axis
    geo = Geodesic(shape.a_axis, flattening)
    lon_0 = phi_0 * 180.0 / math.pi
    lon_1 = phi_1 * 180.0 / math.pi
    lat_0 = 90 - theta_0 * 180.0 / math.pi
    lat_1 = 90 - theta_1 * 180.0 / math.pi
    inv = geo.Inverse(lat_0, lon_0, lat_1, lon_1)
    return inv['s12']
    
    
    
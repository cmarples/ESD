# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 13:50:17 2022

@author: Callum Marples

Functions for computing the taxicab (i.e. rectilinear) distance between two
points on a sphere or spheroid.
"""

import math

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r
def taxicab_distance_sphere(r, theta_0, phi_0, theta_1, phi_1):
    d_theta = r * math.fabs(theta_1-theta_0)
    if math.fabs(theta_0 - math.pi) > math.fabs(theta_1 - math.pi):
        sin_theta = math.sin(theta_0)
    else:
        sin_theta = math.sin(theta_1)
    d_phi = r * sin_theta * math.fabs(phi_1-phi_0)
    return d_theta + d_phi
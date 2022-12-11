# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 13:50:17 2022

@author: Callum Marples

Functions for computing the taxicab (i.e. rectilinear) distance between two
points on a sphere or spheroid.
"""

import math
from scipy.special import ellipeinc

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r.
def taxicab_sphere(r, theta_0, phi_0, theta_1, phi_1):
    pi_by_2 = math.pi/2.0
    d_theta = r * math.fabs(theta_1-theta_0)
    if math.fabs(theta_0 - pi_by_2) > math.fabs(theta_1 - pi_by_2):
        sin_theta = math.sin(theta_0)
    else:
        sin_theta = math.sin(theta_1)
    d_phi = r * sin_theta * math.fabs(phi_1-phi_0)
    return d_theta + d_phi

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a spheroid of axes a and b (with b the distinct axis).
def taxicab_spheroid(a, c, theta_0, phi_0, theta_1, phi_1):
    pi_by_2 = math.pi/2.0
    if c > a: # Prolate
        k2 = 1 - a*a/(c*c)
        d_theta = c * math.fabs((ellipeinc(theta_1-pi_by_2, k2) - ellipeinc(theta_0-pi_by_2, k2)))
    else: # Oblate
        k2 = 1 - c*c/(a*a)
        d_theta = a * math.fabs((ellipeinc(theta_1, k2) - ellipeinc(theta_0, k2)))
    
    if math.fabs(theta_0 - pi_by_2) > math.fabs(theta_1 - pi_by_2):
        sin_theta = math.sin(theta_0)
    else:
        sin_theta = math.sin(theta_1)
    d_phi = a * sin_theta * math.fabs(phi_1-phi_0)
    return d_theta + d_phi




import leod
s = leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 1.0)
s.normalise()
start_point = [50.0, 60.0]
end_point = [90.0, 0.0]
conv = math.pi / 180.0 # Degrees to radians.
d_taxi = taxicab_spheroid(s.a_axis, s.c_axis, start_point[0]*conv, start_point[1]*conv, end_point[0]*conv, end_point[1]*conv)

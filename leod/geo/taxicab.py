# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 16:22:35 2022

@author: Cal
"""

import math
import numpy as np
from scipy.special import ellipeinc
from ..intersection import ellipsoid_plane

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r.
def sphere_tcd(r, start, end, out_flag=False):
    pi_by_2 = math.pi/2.0
    d_theta = r * math.fabs(end[0]-start[0])
    if math.fabs(start[0] - pi_by_2) > math.fabs(end[0] - pi_by_2):
        sin_theta = math.sin(start[0])
    else:
        sin_theta = math.sin(end[0])
    d_phi = r * sin_theta * math.fabs(end[1] - start[1])
    if out_flag == True:
        return [d_theta + d_phi, d_theta, d_phi]
    else:
        return d_theta + d_phi

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a spheroid of axes a and b (with b the distinct axis).
def spheroid_tcd(a, c, start, end, out_flag=False):
    pi_by_2 = math.pi/2.0

    k2 = 1.0 - c*c/(a*a)
    d_theta = a * math.fabs((ellipeinc(end[0], k2) - ellipeinc(start[0], k2)))

    if math.fabs(start[0] - pi_by_2) > math.fabs(end[0] - pi_by_2):
        sin_theta = math.sin(start[0])
    else:
        sin_theta = math.sin(end[0])
    d_phi = a * sin_theta * math.fabs(end[1] - start[1])
    if out_flag == True:
        return [d_theta + d_phi, d_theta, d_phi]
    else:
        return d_theta + d_phi

def triaxial_tcd(shape, start, end, out_flag=False):
    # Constant theta distance
    pi_by_2 = math.pi/2.0
    
    if math.fabs(start[0] - pi_by_2) > math.fabs(end[0] - pi_by_2):
        sin_theta = math.sin(start[0])
        cos_phi = math.cos(end[1])
        sin_phi = math.sin(end[1])
    else:
        sin_theta = math.sin(end[0])
        cos_phi = math.cos(start[1])
        sin_phi = math.sin(start[1])
    a_phi = shape.a_axis*sin_theta
    b_phi = shape.b_axis*sin_theta
    k2 = 1 - b_phi*b_phi/(a_phi*a_phi)
    d_phi = a_phi * math.fabs((ellipeinc(end[1]-pi_by_2, k2) - ellipeinc(start[1]-pi_by_2, k2)))
    
    # Constant phi distance
    sin_th_0 = math.sin(start[0])
    cos_th_0 = math.cos(start[0])
    sin_th_1 = math.sin(end[0])
    cos_th_1 = math.cos(end[0])

    if math.fabs(start[0]-end[0]) < 1.0e-15:
        d_theta = 0.0
    else:
        
        
        r = np.array([shape.a_axis*cos_phi*(sin_th_1 - sin_th_0), shape.b_axis*sin_phi*(sin_th_1 - sin_th_0), shape.c_axis*(cos_th_1 - cos_th_0)])
        s = np.array([0.0, 0.0, 2.0*shape.c_axis])
        normal = np.cross(r, s)
        normal /= np.linalg.norm(normal)
        q = np.array([0.0, 0.0, 0.0])
        centre, A, B, rt, st = ellipsoid_plane(shape, normal, q, True)
        k2 = 1.0 - B*B/(A*A)
        d_theta = A * math.fabs((ellipeinc(end[0], k2) - ellipeinc(start[0], k2)))
        
    if out_flag == True:
        return [d_theta + d_phi, d_theta, d_phi]
    else:
        return d_theta + d_phi


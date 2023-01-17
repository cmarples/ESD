# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 16:22:35 2022

@author: Cal
"""

from math import pi, sin, cos, fabs
import numpy as np
from scipy.special import ellipeinc
from ..intersection import ellipsoid_plane

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r.
def sphere_tcd(r, start, end, out_flag=False, is_radians=False):
    
    start_temp = [0.0, 0.0]
    end_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
    for i in range(2):
        start_temp[i] = start[i] * conv
        end_temp[i] = end[i] * conv
        
    pi_by_2 = pi/2.0
    d_theta = r * fabs(end_temp[0]-start_temp[0])
    if fabs(start_temp[0] - pi_by_2) >fabs(end_temp[0] - pi_by_2):
        sin_theta = sin(start_temp[0])
    else:
        sin_theta = sin(end_temp[0])
    d_phi = r * sin_theta * fabs(end_temp[1] - start_temp[1])
    if out_flag == True:
        return [d_theta + d_phi, d_theta, d_phi]
    else:
        return d_theta + d_phi

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a spheroid of axes a and b (with b the distinct axis).
def spheroid_tcd(a, c, start, end, out_flag=False, is_radians=False):
    
    start_temp = [0.0, 0.0]
    end_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
    for i in range(2):
        start_temp[i] = start[i] * conv
        end_temp[i] = end[i] * conv
        
    pi_by_2 = pi/2.0

    k2 = 1.0 - c*c/(a*a)
    d_theta = a * fabs((ellipeinc(end_temp[0], k2) - ellipeinc(start_temp[0], k2)))

    if fabs(start_temp[0] - pi_by_2) > fabs(end_temp[0] - pi_by_2):
        sin_theta = sin(start_temp[0])
    else:
        sin_theta = sin(end_temp[0])
    d_phi = a * sin_theta * fabs(end_temp[1] - start_temp[1])
    if out_flag == True:
        return [d_theta + d_phi, d_theta, d_phi]
    else:
        return d_theta + d_phi

def triaxial_tcd(shape, start, end, out_flag=False, is_radians=False):
    
    start_temp = [0.0, 0.0]
    end_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
    for i in range(2):
        start_temp[i] = start[i] * conv
        end_temp[i] = end[i] * conv
        
    # Constant theta distance
    pi_by_2 = pi/2.0
    
    if fabs(start_temp[0] - pi_by_2) > fabs(end_temp[0] - pi_by_2):
        sin_theta = sin(start_temp[0])
        cos_phi = cos(end_temp[1])
        sin_phi = sin(end_temp[1])
    else:
        sin_theta = sin(end_temp[0])
        cos_phi = cos(start_temp[1])
        sin_phi = sin(start_temp[1])
    a_phi = shape.a_axis*sin_theta
    b_phi = shape.b_axis*sin_theta
    k2 = 1 - b_phi*b_phi/(a_phi*a_phi)
    d_phi = a_phi * fabs((ellipeinc(end_temp[1]-pi_by_2, k2) - ellipeinc(start_temp[1]-pi_by_2, k2)))
    
    # Constant phi distance
    sin_th_0 = sin(start_temp[0])
    cos_th_0 = cos(start_temp[0])
    sin_th_1 = sin(end_temp[0])
    cos_th_1 = cos(end_temp[0])

    if fabs(start_temp[0]-end_temp[0]) < 1.0e-15:
        d_theta = 0.0
    else:
        
        
        r = np.array([shape.a_axis*cos_phi*(sin_th_1 - sin_th_0), shape.b_axis*sin_phi*(sin_th_1 - sin_th_0), shape.c_axis*(cos_th_1 - cos_th_0)])
        s = np.array([0.0, 0.0, 2.0*shape.c_axis])
        normal = np.cross(r, s)
        normal /= np.linalg.norm(normal)
        q = np.array([0.0, 0.0, 0.0])
        centre, A, B, rt, st = ellipsoid_plane(shape, normal, q, True)
        k2 = 1.0 - B*B/(A*A)
        d_theta = A * fabs((ellipeinc(end_temp[0], k2) - ellipeinc(start_temp[0], k2)))
        
    if out_flag == True:
        return [d_theta + d_phi, d_theta, d_phi]
    else:
        return d_theta + d_phi


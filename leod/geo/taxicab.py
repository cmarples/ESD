# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 16:22:35 2022

@author: Cal
"""

import math
import numpy as np
from scipy.special import ellipeinc

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r.
def sphere_td(r, start, end):
    pi_by_2 = math.pi/2.0
    d_theta = r * math.fabs(end[0]-start[0])
    if math.fabs(start[0] - pi_by_2) > math.fabs(end[0] - pi_by_2):
        sin_theta = math.sin(start[0])
    else:
        sin_theta = math.sin(end[0])
    d_phi = r * sin_theta * math.fabs(end[1] - start[1])
    return d_theta + d_phi

# Calculate taxicab distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a spheroid of axes a and b (with b the distinct axis).
def spheroid_td(a, c, start, end):
    pi_by_2 = math.pi/2.0
    if c > a: # Prolate
        k2 = 1.0 - a*a/(c*c)
        d_theta = c * math.fabs((ellipeinc(end[0]-pi_by_2, k2) - ellipeinc(start[0]-pi_by_2, k2)))
    else: # Oblate
        k2 = 1.0 - c*c/(a*a)
        d_theta = a * math.fabs((ellipeinc(end[0], k2) - ellipeinc(start[0], k2)))
    
    if math.fabs(start[0] - pi_by_2) > math.fabs(end[0] - pi_by_2):
        sin_theta = math.sin(start[0])
    else:
        sin_theta = math.sin(end[0])
    d_phi = a * sin_theta * math.fabs(end[1] - start[1])
    return d_theta + d_phi


def triaxial_td(a, b, c, start, end):
    pi_by_2 = math.pi/2.0
    
    sin_th_0 = math.sin(start[0])
    cos_th_0 = math.cos(start[0])
    sin_th_1 = math.sin(end[0])
    cos_th_1 = math.cos(end[0])
    
    # Constant theta distance
    if math.fabs(start[0] - pi_by_2) > math.fabs(end[0] - pi_by_2):
        sin_theta = math.sin(start[0])
    else:
        sin_theta = math.sin(end[0])
    a_phi = a*sin_theta
    b_phi = b*sin_theta
    k2 = 1 - b_phi*b_phi/(a_phi*a_phi)
    d_phi = a_phi * math.fabs((ellipeinc(start[1], k2) - ellipeinc(start[0], k2)))

    # Constant phi distance
    x0 = a_phi * math.cos(start[1])
    y0 = a_phi * math.sin(start[1])
    x1 = a_phi * math.cos(end[1])
    y1 = a_phi * math.sin(end[1])
    if math.fabs(math.sqrt(x0*x0+y0*y0)) > math.fabs(math.sqrt(x1*x1+y1*y1)):
        phi = start[1]
        cos_phi = math.cos(start[1])
        sin_phi = math.sin(start[1])
    else:
        phi = end[1]
        cos_phi = math.cos(end[1])
        sin_phi = math.sin(end[1])
    
    r = np.array([a*cos_phi*(sin_th_1 - sin_th_0), b*sin_phi*(sin_th_1 - sin_th_0), c*(cos_th_1 - cos_th_0)])
    s = np.array([0.0, 0.0, 2.0*c])
    Dr = np.array([r[0]/a, r[1]/b, r[2]/c])
    Ds = np.array([0.0, 0.0, 2.0])
    w = 0.5*math.atan2(2.0*np.dot(Dr, Ds), np.dot(Dr, Dr) - np.dot(Ds, Ds))
    cw = math.cos(w)
    sw = math.sin(w)
    r_tilde = cw*r + sw*s
    s_tilde = -sw*r + cw*s
    Dr = np.array([r_tilde[0]/a, r_tilde[1]/b, r_tilde[2]/c])
    Ds = np.array([s_tilde[0]/a, s_tilde[1]/b, s_tilde[2]/c])
    a_theta = math.sqrt(1.0 / np.dot(Dr, Dr))
    b_theta = math.sqrt(1.0 / np.dot(Ds, Ds))
    
    k2_th = 1 - b_theta*b_theta/(a_theta*a_theta)
    d_theta = a_theta * math.fabs((ellipeinc(end[1], k2_th) - ellipeinc(start[1], k2_th)))
    
    return d_theta + d_phi
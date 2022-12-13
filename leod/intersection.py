# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:09:40 2022

@author: Cal
"""

import math
import numpy as np

def ellipsoid_plane(shape, n_vec, q, out_flag=False):
    
    Dq = np.array([q[0]/shape.a_axis, q[1]/shape.b_axis, q[2]/shape.c_axis])
    
    # Find in-plane vectors, r and s, such that the two are mutually orthonormal.
    if math.fabs(n_vec[0]) < 1.0e-15 and math.fabs(n_vec[1]) < 1.0e-15:
        r = np.array([1.0, 0.0, 0.0])
        s = np.array([0.0, 1.0, 0.0])
    else:
        r = np.array([n_vec[1], -n_vec[0], 0.0])
        s = np.cross(n_vec, r)
        # Normalise to one
        r /= np.linalg.norm(r)
        s /= np.linalg.norm(s)
    
    
    Dr = np.array([r[0]/shape.a_axis, r[1]/shape.b_axis, r[2]/shape.c_axis])
    Ds = np.array([s[0]/shape.a_axis, s[1]/shape.b_axis, s[2]/shape.c_axis])
    
    if np.dot(Dr, Ds) != 0.0:
        # Rotate by w to make Dr.Ds=0
        w = 0.5*math.atan2(2.0*np.dot(Dr, Ds), np.dot(Dr, Dr) - np.dot(Ds, Ds))
        cw = math.cos(w)
        sw = math.sin(w)
        r_tilde = cw*r + sw*s
        s_tilde = -sw*r + cw*s
        Dr = np.array([r_tilde[0]/shape.a_axis, r_tilde[1]/shape.b_axis, r_tilde[2]/shape.c_axis])
        Ds = np.array([s_tilde[0]/shape.a_axis, s_tilde[1]/shape.b_axis, s_tilde[2]/shape.c_axis])
    
    # Dot products
    DrDr = np.dot(Dr, Dr)
    DsDs = np.dot(Ds, Ds)
    DqDr = np.dot(Dq, Dr)
    DqDs = np.dot(Dq, Ds)
    DqDq = np.dot(Dq, Dq)
    
    # Ellipse centre and axes
    center = [-DqDr/DrDr, -DqDs/DsDs]
    d = DqDq - DqDr*DqDr/(DrDr) - DqDs*DqDs/(DsDs)
    a_ellipse = math.sqrt((1.0-d) / DrDr)
    b_ellipse = math.sqrt((1.0-d) / DsDs)   
    
    if out_flag == True:
        return [center, a_ellipse, b_ellipse, r, s]
    else:
        return [center, a_ellipse, b_ellipse]
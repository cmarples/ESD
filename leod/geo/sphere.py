# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:25:04 2022

@author: Callum Marples
"""

import math
import numpy as np

# Great circle distance.
# Calculate shortest distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r
def gc_dist(r, start, end):
    """! Shortest distance on a sphere of radius \f$r\f$, using a great circle arc length.
    @param r : float \n
        The sphere radius.
    @param start : list of floats \n
        The start point in polar coordinates, \f$[\theta_0, \phi_0]\f$.
    @param end : list of floats \n
        The end point in polar coordinates, \f$[\theta_1, \phi_1]\f$.
    @return The shortest distance between start and end.
    """
    return r * math.acos( math.cos(start[0])*math.cos(end[0]) +
                          math.sin(start[0])*math.sin(end[0]) *
                          math.cos(end[1] - start[1]) )

# Get shortest path between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r
def gc_path(r, start, end):
    """! Shortest path on a sphere of radius \f$r\f$, using a great circle arc length.
    @param r : float \n
        The sphere radius.
    @param start : list of floats \n
        The start point in polar coordinates, \f$[\theta_0, \phi_0]\f$.
    @param end : list of floats \n
        The end point in polar coordinates, \f$[\theta_1, \phi_1]\f$.
    @return The shortest path between start and end.
    """
    theta_0 = start[0]
    phi_0 = start[1]
    theta_1 = end[0]
    phi_1 = end[1]
    phi_vals = np.linspace(phi_0, phi_1, 100)
    if abs(phi_1 - phi_0) > 1e-15:
        T = math.tan(theta_1) / math.tan(theta_0)
        phi_c = math.atan2(math.cos(phi_0) - T*math.cos(phi_1), T*math.sin(phi_1) - math.sin(phi_0))
        a = 1 / (math.tan(theta_0) * math.cos(phi_0 - phi_c))
        cot_th_vals = a * np.cos(phi_vals - phi_c)
        th_vals = np.arctan2(1, cot_th_vals)
    else:
        th_vals = np.linspace(theta_0, theta_1, 100)
    
    path_positions = np.zeros((len(th_vals), 3)) 
    # Cartesian coordinates
    x = r * np.sin(th_vals) * np.cos(phi_vals)
    y = r * np.sin(th_vals) * np.sin(phi_vals)
    z = r * np.cos(th_vals)
    for i in range(len(th_vals)):
        path_positions[i] = [x[i], y[i], z[i]]
    return path_positions
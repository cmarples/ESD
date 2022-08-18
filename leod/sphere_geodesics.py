# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:25:04 2022

@author: Callum Marples
"""

import math
import numpy as np

# Calculate shortest distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r
def great_circle_distance(r, theta_0, phi_0, theta_1, phi_1):
    return r * math.acos( math.cos(theta_0)*math.cos(theta_1) +
                          math.sin(theta_0)*math.sin(theta_1) *
                          math.cos(phi_0 - phi_1) )

# Get shortest path between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r
def great_circle_path(r, theta_0, phi_0, theta_1, phi_1):
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
"""! 
@brief Spherical geodesic (great circle) routines.
@file sphere.py
@author Callum Marples
- Created by Callum Marples on 16/08/2022.
- Last modified on 31/01/2022.
"""

from math import pi, sin, cos, tan, acos, atan2
import numpy as np

# Great circle distance.
# Calculate shortest distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r
def gc_dist(r, start, end, is_radians=False):
    """! Shortest distance on a sphere of radius \f$r\f$, using a great circle arc length.
    @param r : float \n
        The sphere radius.
    @param start : list of floats \n
        The start point in polar coordinates, \f$[\theta_0, \phi_0]\f$.
    @param end : list of floats \n
        The end point in polar coordinates, \f$[\theta_1, \phi_1]\f$.
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return The shortest distance between start and end.
    """
    
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
        
    return r * acos( cos(start_temp[0])*cos(end_temp[0]) +
                     sin(start_temp[0])*sin(end_temp[0]) *
                     cos(end_temp[1] - start_temp[1]) )

# Get shortest path between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a sphere of radius r
def gc_path(r, start, end, is_radians=False):
    """! Shortest path on a sphere of radius \f$r\f$, using a great circle arc length.
    @param r : float \n
        The sphere radius.
    @param start : list of floats \n
        The start point in polar coordinates, \f$[\theta_0, \phi_0]\f$.
    @param end : list of floats \n
        The end point in polar coordinates, \f$[\theta_1, \phi_1]\f$.
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return The shortest path between start and end.
    """
    
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
        
    theta_0 = start_temp[0]
    phi_0 = start_temp[1]
    theta_1 = end_temp[0]
    phi_1 = end_temp[1]
    
    phi_vals = np.linspace(phi_0, phi_1, 100)
    if abs(phi_1 - phi_0) > 1e-15:
        T = tan(theta_1) / tan(theta_0)
        phi_c = atan2(cos(phi_0) - T*cos(phi_1), T*sin(phi_1) - sin(phi_0))
        a = 1 / (tan(theta_0) * cos(phi_0 - phi_c))
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
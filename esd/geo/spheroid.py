"""! 
@brief Spheroid geodesic routine.
@file spheroid.py
@author Callum Marples
- Created by Callum Marples on 5/09/2022.
- Last modified on 31/01/2022.
"""

from math import pi
from geographiclib.geodesic import Geodesic

# Call GeographicLib
# Calculate shortest distance between (theta_0, phi_0) and (theta_1, phi_1)
# on the surface of a spheroid of axis lengths (a, a, c), as given in an
# EllipsoidShape object
def geo_dist(shape, start, end, is_radians=False):
    """! Calculate spheroid geodesic distance, using GeographicLib \cite Karney2019.
    To use this routine, GeographicLib must be installed!
    @param shape : EllipsoidShape \n
        The spheroid (must input an EllipsoidShape object such that a_axis and b_axis are equal).
    @param start : list of floats \n
        The start point in polar coordinates, \f$[\theta_0, \phi_0]\f$.
    @param end : list of floats \n
        The end point in polar coordinates, \f$[\theta_1, \phi_1]\f$.
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordinate conversion to
        radians is performed. Defaults to False.
    @return The geodesic distance between start and end.
    """
    start_temp = [0.0, 0.0]
    end_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = pi / 180.0
        conv_inv = 1.0 / conv
    else:
        conv = 1.0
        conv_inv = 1.0
    for i in range(2):
        start_temp[i] = start[i] * conv
        end_temp[i] = end[i] * conv
        
    theta_0 = start_temp[0]
    phi_0 = start_temp[1]
    theta_1 = end_temp[0]
    phi_1 = end_temp[1]
    
    flattening = (shape.a_axis - shape.c_axis) / shape.a_axis
    geo = Geodesic(shape.a_axis, flattening)
    lon_0 = phi_0 * conv_inv
    lon_1 = phi_1 * conv_inv
    lat_0 = 90 - theta_0 * conv_inv
    lat_1 = 90 - theta_1 * conv_inv
    inv = geo.Inverse(lat_0, lon_0, lat_1, lon_1)
    return inv['s12']
    
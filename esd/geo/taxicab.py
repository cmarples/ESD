"""! 
@brief Taxicab (rectilinear) distance routines.
@file taxicab.py
@author Callum Marples
- Created by Callum Marples on 12/12/2022.
- Last modified on 31/01/2022.
"""

from math import pi, sin, cos, fabs
import numpy as np
from scipy.special import ellipeinc
from ..intersection import ellipsoid_plane

def sphere_tcd(r, start, end, out_flag=False, is_radians=False):
    """! Calculate shortest taxicab distance between two points on the surface of a sphere.
    The path used is equivalent to that described by Bayar and Kaya \cite Bayar2005.
    @param r : float \n
        Radius of the sphere.
    @param start : list of floats (2 elements) \n
        The start point, in \f$(\theta, \phi)\f$ coordinates.
    @param end : list of floats (2 elements) \n
        The end point, in \f$(\theta, \phi)\f$ coordinates.
    @param out_flag : bool (optional)
        If True, the individual \f$\theta\f$ and \f$\phi\f$ distances are 
        included in the output. Defaults to False. 
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return float or 3 element list of float \n
        If out_flag False, outputs the taxicab distance from start to end.
        If out_flag True, outputs a list containing taxicab distance, 
        constant \f$\theta\f$ distance and constant \f$\phi\f$ distance respectively.
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



def spheroid_tcd(a, c, start, end, out_flag=False, is_radians=False):
    """! Calculate shortest taxicab distance between two points on the surface of a spheroid.
    The input spheroid may be oblate (\f$a>c\f$) or prolate (\f$a<c\f$).
    The constant \f$\theta\f$ distance is computed using elliptic integrals.
    @param a : float \n
        The repeated axis, \f$a=b\f$, of the spheroid.
    @param c : float \n
        The distinct axis, \f$c\f$, of the spheroid.
    @param start : list of floats (2 elements) \n
        The start point, in \f$(\theta, \phi)\f$ coordinates.
    @param end : list of floats (2 elements) \n
        The end point, in \f$(\theta, \phi)\f$ coordinates.
    @param out_flag : bool (optional)
        If True, the individual \f$\theta\f$ and \f$\phi\f$ distances are 
        included in the output. Defaults to False. 
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return float or 3 element list of float \n
        If out_flag False, outputs the taxicab distance from start to end.
        If out_flag True, outputs a list containing taxicab distance, 
        constant \f$\theta\f$ distance and constant \f$\phi\f$ distance respectively.
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
    """! Calculate shortest taxicab distance between two points on the surface of a triaxial ellipsoid.
    Both the constant \f$\theta\f$ and \f$\phi\f$ distances are computed using elliptic integrals.
    The appropriate ellipse for the constant \f$\theta\f$ distance is found by 
    finding the intersection between the ellipsoid and a plane containing the start point, 
    end point and the ellipsoid centre.
    @param shape : EllipsoidShape \n
        The ellipsoid.
    @param start : list of floats (2 elements) \n
        The start point, in \f$(\theta, \phi)\f$ coordinates.
    @param end : list of floats (2 elements) \n
        The end point, in \f$(\theta, \phi)\f$ coordinates.
    @param out_flag : bool (optional)
        If True, the individual \f$\theta\f$ and \f$\phi\f$ distances are 
        included in the output. Defaults to False. 
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return float or 3 element list of float \n
        If out_flag False, outputs the taxicab distance from start to end.
        If out_flag True, outputs a list containing taxicab distance, 
        constant \f$\theta\f$ distance and constant \f$\phi\f$ distance respectively.
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
        # Determine the intersecting ellipse for a plane containing the start 
        # point, end point and ellipsoid centre.
        r = np.array([shape.a_axis*cos_phi*(sin_th_1 - sin_th_0), shape.b_axis*sin_phi*(sin_th_1 - sin_th_0), shape.c_axis*(cos_th_1 - cos_th_0)])
        s = np.array([0.0, 0.0, 2.0*shape.c_axis])
        normal = np.cross(r, s)
        normal /= np.linalg.norm(normal)
        q = np.array([0.0, 0.0, 0.0])
        A, B = ellipsoid_plane(shape, normal, q)
        # Using this ellipse, determine constant \f$\theta\f$ distance.
        k2 = 1.0 - B*B/(A*A)
        d_theta = A * fabs((ellipeinc(end_temp[0], k2) - ellipeinc(start_temp[0], k2)))
        
    if out_flag == True:
        return [d_theta + d_phi, d_theta, d_phi]
    else:
        return d_theta + d_phi


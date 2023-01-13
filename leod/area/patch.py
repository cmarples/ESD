# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 12:25:25 2023

@author: Callum Marples

Routines for patch areas of a spheroid/ellipsoid.
"""

from math import sqrt, pi, fabs, cos, asin, asinh
from scipy.special import ellipeinc
from scipy.integrate import romberg
from numpy import array

### Analytical
def sphere(r, th_0, th_1, ph_0, ph_1):
    """! Surface area of a spherical patch; between polar angles \f$\theta_0\f$ and \f$\theta_1\f$,
         and azimuthal angles \f$\phi_0\f$ and \f$\phi_1\f$.
         The analytical expression for spherical patch area is used.
    @param r : float \n
        The radius of the sphere.
    @param th_0 : float \n
        Start angle, \f$\theta_0\f$ of the patch.
    @param th_1 : float \n
        End angle, \f$\theta_1\f$ of the patch.
    @param ph_0 : float \n
        Start angle, \f$\theta_0\f$ of the patch.
    @param ph_1 : float \n
        End angle, \f$\theta_1\f$ of the patch.
    @return double \n
        Surface area of the spherical patch.
    """
    return r*r*(ph_1 - ph_0)*(cos(th_0) - cos(th_1))

def spheroid_band(shape, th_0, th_1):
    """! Surface area of a spheroidal band, between polar angles \f$\theta_0\f$ and \f$\theta_1\f$.
         The shape is assumed to have repeated axis \f$a\f$, and distinct axis \f$c\f$ (i.e. \f$a=b\f$).
         In this routine, the analytical expressions for spheroid band areas are used.
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @param th_0 : float \n
        Start angle, \f$\theta_0\f$ of the band.
    @param th_1 : float \n
        End angle, \f$\theta_1\f$ of the band.
    @return double \n
        Surface area of the spheroidal band.
    """
    if th_0 < 0.0 or th_0 > pi:
        raise ValueError(f"th_0 not in range [0, pi]")
    if th_1 < 0.0 or th_1 > pi:
        raise ValueError(f"th_1 not in range [0, pi]")
    if shape.is_spheroid() == False:
        raise ValueError(f"Input shape is not a spheroid.")
    q = 1.0 - shape.a_axis*shape.a_axis/(shape.c_axis*shape.c_axis)
    # z coordinates of the band boundaries.
    cos_0 = cos(th_0)
    cos_1 = cos(th_1)
    if shape.a_axis < shape.c_axis: # Prolate
        q = sqrt(q)
        return pi*shape.a_axis*shape.c_axis*( (asin(q*cos_0) - asin(q*cos_1))/q + 
                                              cos_0*sqrt(1.0 - q*cos_0*cos_0) -
                                              cos_1*sqrt(1.0 - q*cos_1*cos_1) )
    else:                           # Oblate
        q = sqrt(-q)
        return pi*shape.a_axis*shape.c_axis*( (asinh(q*cos_0) - asinh(q*cos_1))/q + 
                                              cos_0*sqrt(1.0 - q*cos_0*cos_0) -
                                              cos_1*sqrt(1.0 - q*cos_1*cos_1) )
    
def spheroid(shape, th_0, th_1, ph_0, ph_1):
    """! Surface area of a spheroidal patch; between polar angles \f$\theta_0\f$ and \f$\theta_1\f$,
         and azimuthal angles \f$\phi_0\f$ and \f$\phi_1\f$.
         The shape is assumed to have repeated axis \f$a\f$, and distinct axis \f$c\f$ (i.e. \f$a=b\f$).
         In this routine, the analytical expressions for spheroid band areas are used.
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @param th_0 : float \n
        Start angle, \f$\theta_0\f$ of the patch.
    @param th_1 : float \n
        End angle, \f$\theta_1\f$ of the patch.
    @param ph_0 : float \n
        Start angle, \f$\theta_0\f$ of the patch.
    @param ph_1 : float \n
        End angle, \f$\theta_1\f$ of the patch.
    @return double \n
        Surface area of the spheroidal patch.
    """
    if ph_0 < 0.0 or ph_0 > 2.0*pi:
        raise ValueError(f"ph_0 not in range [0, 2*pi]")
    if ph_1 < 0.0 or ph_1 > 2.0*pi:
        raise ValueError(f"ph_1 not in range [0, 2*pi]")
    band_area = spheroid_band(shape, th_0, th_1)
    delta_ph = fabs(ph_1 - ph_0)
    # Account for periodicity of phi
    if delta_ph > pi:
        if ph_0 > pi:
            ph_0 -= 2.0*pi
        else:
            ph_1 -= 2.0*pi
        delta_ph = fabs(ph_1 - ph_0)
    return band_area * delta_ph/(2.0*pi)
    
    
def ellipsoid(shape, th_0, th_1, ph_0, ph_1):
    """! Surface area of an ellipsoidal patch; between polar angles \f$\theta_0\f$ and \f$\theta_1\f$,
         and azimuthal angles \f$\phi_0\f$ and \f$\phi_1\f$.
         Romberg integration is used to numerically obtain the area.
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @param th_0 : float \n
        Start angle, \f$\theta_0\f$ of the patch.
    @param th_1 : float \n
        End angle, \f$\theta_1\f$ of the patch.
    @param ph_0 : float \n
        Start angle, \f$\theta_0\f$ of the patch.
    @param ph_1 : float \n
        End angle, \f$\theta_1\f$ of the patch.
    @return double \n
        Surface area of the ellipsoidal patch.
    """    
    if shape.is_sphere() == True:
        return sphere(shape.a_axis, th_0, th_1, ph_0, ph_1)
    else:
        t0 = cos(th_0)
        t1 = cos(th_1)
        if shape.a_axis > shape.b_axis:
            b2 = shape.b_axis*shape.b_axis
            g_frac = b2/(shape.c_axis*shape.c_axis) - 1.0
            m_frac = 1.0 - b2/(shape.a_axis*shape.a_axis)
            pre_fac = shape.a_axis*shape.c_axis
            ph_0 -= 0.5*pi
            ph_1 -= 0.5*pi
        else:
            a2 = shape.a_axis*shape.a_axis
            g_frac = a2/(shape.c_axis*shape.c_axis) - 1.0
            m_frac = 1.0 - a2/(shape.b_axis*shape.b_axis)
            pre_fac = shape.b_axis*shape.c_axis
        def integrand(x):
            g = 1.0 + g_frac*x*x
            m = m_frac*(1.0 - x*x) / g
            return sqrt(g)*(ellipeinc(ph_1, m) - ellipeinc(ph_0, m))
        return pre_fac*romberg(integrand, t1, t0)
    
    
def grid(shape, grid):
    """! Numerically calculate surface areas of all patches in a grid.
         This routine uses the 'ellipsoid' routine to calculate areas.
         The 8-fold symmetry of the triaxial ellipsoid is exploited to reduce
         the number of 'ellipsoid' calls.
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @param grid : Grid \n
        The grid of vertices/patches on the ellipsoid surface.
    @return NumPy array of floats \n
        Surface areas of each ellipsoidal patch in the input grid.
    """    
    # Initialise output array
    patch_areas = [0.0] * grid.no_vertices
    patch_areas = array(patch_areas) 
    # Determine symmetry of the grid (i.e. no_theta odd and/or no_phi multiple of 4)
    theta_end = int((grid.no_theta - 1) / 2)
    phi_end = int((grid.no_phi - 2) / 4)
    if grid.no_theta % 2 == 0:
        is_theta_odd = False
    else:
        is_theta_odd = True
    m = 0
    if grid.no_phi % 4 == 0:
        is_phi_4 = True
        m = (grid.no_phi - 4) / 4 # Number of patches between a and b axes.
        phi_end = int(m) + 1
    else:
        is_phi_4 = False
        
    theta = 0.0 # theta value at current patch
    k = 0 # Scalar index of current patch
    for i in range(theta_end + 1):   # theta loop
        phi = 0.0 # phi value at current patch
        for j in range(phi_end + 1): # phi loop
            # Calculate patch area
            if k == 0: # Polar cap areas
                patch_areas[0] = 4.0*ellipsoid(shape, 0.0, 0.5*grid.delta_theta, 0.0, 0.5*pi)
                patch_areas[grid.no_vertices - 1] = patch_areas[0]
                k += 1
                break # Exit phi loop and proceed to next theta
            else:
                lmda = 1.0
                # theta limits
                if (i == theta_end and is_theta_odd==True):
                    th_0 = 0.5*(pi - grid.delta_theta)
                    th_1 = 0.5*pi
                    lmda *= 2.0
                else:
                    th_0 = theta - 0.5*grid.delta_theta
                    th_1 = theta + 0.5*grid.delta_theta
                # phi limits
                if (j == phi_end and is_phi_4==True):
                    ph_0 = 0.5*(pi - grid.delta_phi)
                    ph_1 = 0.5*pi
                    lmda *= 2.0
                elif j == 0: ##### NEW!
                    ph_0 = 0.0
                    ph_1 = 0.5*grid.delta_phi
                    lmda *= 2.0
                else:
                    ph_0 = phi - 0.5*grid.delta_phi
                    ph_1 = phi + 0.5*grid.delta_phi
                    
                # Perform patch area calculation.
                patch_areas[k] = lmda * ellipsoid(shape, th_0, th_1, ph_0, ph_1)
                
                # Use symmetry to obtain areas of equivalent patches.
                if is_phi_4 == False:
                    if j == 0:
                        patch_areas[grid.get_vertex_index(i, grid.no_phi/2)] = patch_areas[k]
                        if (is_theta_odd==False or i != theta_end):
                            patch_areas[grid.get_vertex_index(grid.no_theta - 1 - i, 0)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta - 1 - i, int(grid.no_phi/2))] = patch_areas[k]
                    else:
                        patch_areas[grid.get_vertex_index(i, phi_end + j)] = patch_areas[k]
                        patch_areas[grid.get_vertex_index(i, int(grid.no_phi/2) + j)] = patch_areas[k]
                        patch_areas[grid.get_vertex_index(i, int(grid.no_phi/2) + phi_end + j)] = patch_areas[k]
                        if (is_theta_odd == False or i != theta_end):
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, phi_end + j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, int(grid.no_phi/2) + j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, int(grid.no_phi/2) + phi_end + j)] = patch_areas[k]
                else:
                    if j == 0:
                        patch_areas[grid.get_vertex_index(i, int(grid.no_phi/2))] = patch_areas[k]
                        if (is_theta_odd==False or i != theta_end):
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, int(grid.no_phi/2) + j)] = patch_areas[k]
                    elif j == phi_end: ### ADDED!
                        patch_areas[grid.get_vertex_index(i, int(grid.no_phi/2) + j)] = patch_areas[k]
                        if (is_theta_odd==False or i != theta_end):
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, int(grid.no_phi/2) + j)] = patch_areas[k]
                    
                    else:
                        patch_areas[grid.get_vertex_index(i, 2*phi_end-j)] = patch_areas[k]
                        patch_areas[grid.get_vertex_index(i, 2*phi_end+j)] = patch_areas[k]
                        patch_areas[grid.get_vertex_index(i, 4*phi_end-j)] = patch_areas[k]
                        if (is_theta_odd==False or i != theta_end):
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, 2*phi_end-j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, 2*phi_end+j)] = patch_areas[k]
                            patch_areas[grid.get_vertex_index(grid.no_theta-1-i, 4*phi_end-j)] = patch_areas[k]
                # Advance counter k (patch index)
                if j == phi_end:
                    k = grid.get_vertex_index(i+1, 0)
                else:
                    k += 1
            phi += grid.delta_phi
        theta += grid.delta_theta
    return patch_areas
                                
 
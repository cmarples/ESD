# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 12:19:46 2023

@author: Callum Marples

Routines for the full surface area of a spheroid/ellipsoid.
"""

from math import sqrt, pi, asin, acos, asinh
from scipy.special import ellipeinc, ellipkinc, ellipe
from scipy.integrate import romberg

### Analytical
def sphere(r):
    """! Surface area of a sphere of radius \f$r\f$.
    @param r : float \n
        Radius of the sphere.
    @return double \n
        Surface area, \f$4\pi r^2\f$.
    """
    return 4*pi*r*r

def spheroid(shape):
    """! Surface area of a spheroid. The shape is assumed to have repeated axis \f$a\f$, and distinct axis \f$c\f$ (i.e. \f$a=b\f$).
         In this routine, the analytical expressions for spheroid area are used.
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @return double \n
        Surface area of the spheroid.
    """
    # Check that the shape really is a spheroid
    if shape.is_spheroid() == False:
        raise ValueError(f"Input shape is not a spheroid.")
    q = 1.0 - shape.a_axis*shape.a_axis/(shape.c_axis*shape.c_axis)
    if shape.a_axis < shape.c_axis: # Prolate
        q = sqrt(q)
        return 2.0*pi*shape.a_axis*( shape.a_axis + shape.c_axis*asin(q)/q )
    else:                           # Oblate
        q = sqrt(-q)
        return 2.0*pi*shape.a_axis*( shape.a_axis + shape.c_axis*asinh(q)/q )
    
def thomsen(shape, p=1.6075):
    """! Approximate surface area of an ellipsoid using Thomsen's approximation.
         See http://www.numericana.com/answer/ellipsoid.htm#thomsen
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @param p : float (optional) \n
        Parameter used in the approximate formula, defaults to 1.6075.
    @return double \n
        Approximate surface area of the spheroid.
    """
    x = ((shape.a_axis*shape.b_axis)**p + (shape.a_axis*shape.c_axis)**p + (shape.b_axis*shape.c_axis)**p) / 3.0
    return 4.0*pi*x**(1.0/p)

### Numerical
def legendre(shape):
    """! Surface area of an ellipsoid using Legendre's formula.
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @return double \n
        Surface area of the spheroid.
    """
    cos_2_phi = shape.c_axis / shape.a_axis # Here stores cos(phi)
    phi = acos(cos_2_phi)
    cos_2_phi *= cos_2_phi # Now stores cos^2(phi)
    sin_2_phi = 1.0 - cos_2_phi
    a2 = shape.a_axis*shape.a_axis
    b2 = shape.b_axis*shape.b_axis
    c2 = shape.c_axis*shape.c_axis
    m = (a2*(b2-c2)) / (b2*(a2-c2))
    return 2.0*pi*( c2 + shape.a_axis*shape.b_axis*(sin_2_phi*ellipeinc(phi, m) + cos_2_phi*ellipkinc(phi, m)) / sqrt(sin_2_phi) )
    

def numerical(shape):
    """! Surface area of an ellipsoid using Romberg integration.
    @param shape : EllipsoidShape \n
        The ellipsoid shape object.
    @return double \n
        Surface area of the spheroid.
    """
    delta = 1.0 - shape.c_axis*shape.c_axis/(shape.a_axis*shape.a_axis)
    epsln = 1.0 - shape.c_axis*shape.c_axis/(shape.b_axis*shape.b_axis)
    # Define integrand
    def integrand(x):
        x2 = x*x
        m = epsln*(1.0-x2) / (1.0 - delta*x2)
        return sqrt(1.0 - delta*x2) * ellipe(m)
    return 8.0*shape.a_axis*shape.b_axis*romberg(integrand, 0.0, 1.0)

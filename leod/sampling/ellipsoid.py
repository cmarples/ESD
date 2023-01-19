# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:59:56 2023

@author: Callum Marples

Random sampling of the surface of a sphere.
"""

from math import sqrt, pi, cos, sin, exp, fabs
from numpy import array
from numpy.linalg import norm
from numpy.random import Generator, PCG64, MT19937
from .sphere import Marsaglia, Trig

class EllipsoidBase:
    
    def initialise(self, shape, rng_type="MT", seed=12345):
        """! Base initialisation for all ellipsoid sampling classes.
        @param shape : EllipsoidShape \n
            The ellipsoid whose surface is to be sampled.
        @param rng_type : str (optional) \n
            Specifies the random number generator to use. There are two options: \n
            - "MT"  : Mersenne twister.
            - "PCG" : PCG64.
            The default setting is "MT".
        @param seed : int (optional) \n
            The seed of the random number generator, defaults to 12345.
        """
        if rng_type == "PCG":
            self.rng = Generator(PCG64(seed))
        else:
            self.rng = Generator(MT19937(seed))
        self.shape = shape
        self.no_attempts = 0


class NaiveScale(EllipsoidBase):
    
    def __init__(self, shape, rng_type="MT", seed=12345):
        """! Initialiser for the naive scaling ellipsoid sampler.
        @param shape : EllipsoidShape \n
            The ellipsoid whose surface is to be sampled.
        @param rng_type : str (optional) \n
            Specifies the random number generator to use. There are two options: \n
            - "MT"  : Mersenne twister (default setting).
            - "PCG" : PCG64.
        @param seed : int (optional) \n
            The seed of the random number generator, defaults to 12345.
        """
        self.initialise(shape, rng_type, seed)
        self.sphere_sampler = Marsaglia(1.0)
        self.sphere_sampler.set_rng(self.rng) # Make both sphere and ellipsoid samplers use the same rng object.
        
    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of an ellipsoid, 
        using the naive scaling method.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """        
        p = self.sphere_sampler.random_surface_point() # Generate sphere point.
        p[0] *= self.shape.a_axis
        p[1] *= self.shape.b_axis # Anisotropic scaling
        p[2] *= self.shape.c_axis
        self.no_attempts += 1
        return p

        
class GradRej(EllipsoidBase):
    
    def __init__(self, shape, rng_type="MT", seed=12345, sphere_type="Trig"):
        """! Initialiser for the gradient rejection ellipsoid sampler.
        @param shape : EllipsoidShape \n
            The ellipsoid whose surface is to be sampled.
        @param rng_type : str (optional) \n
            Specifies the random number generator to use. There are two options: \n
            - "MT"  : Mersenne twister (default setting).
            - "PCG" : PCG64.
        @param seed : int (optional) \n
            The seed of the random number generator, defaults to 12345.
        @param sphere_type : str \n
            Specifies the sphere sampler to use
            - "Trig"      : Trigonometric method (default setting).
            - "Marsaglia" : Marsaglia's method.
        """
        self.initialise(shape, rng_type, seed)
        if sphere_type == "Marsaglia":
            self.sphere_sampler = Marsaglia(1.0)
        else:
            self.sphere_sampler = Trig(1.0)
        self.sphere_sampler.set_rng(self.rng) # Make both sphere and ellipsoid samplers use the same rng object.
        # Pre-calculate frequently used values
        self.g_min = shape.c_axis
        if shape.b_axis < self.g_min:
            self.g_min = shape.b_axis
        if shape.a_axis < self.g_min:
            self.g_min = shape.a_axis
        self.a4 = 1.0 / shape.a_axis ** 4
        self.b4 = 1.0 / shape.b_axis ** 4
        self.c4 = 1.0 / shape.c_axis ** 4
        
    def random_surface_point(self, out_type="carts"):
        """! Generates a uniformly random point on the surface of an ellipsoid, 
        using the gradient rejection method.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """        
        accept = False
        while accept == False:
            p = self.sphere_sampler.random_surface_point() # Generate sphere point.
            p[0] *= self.shape.a_axis
            p[1] *= self.shape.b_axis # Anisotropic scaling
            p[2] *= self.shape.c_axis
            g = self.g_min * sqrt(self.a4*p[0]*p[0] + self.b4*p[1]*p[1] + self.c4*p[2]*p[2]) # Acceptance probability
            u = self.rng.random() # 'Roll dice'
            if u <= g:
                accept = True
            self.no_attempts += 1
        
        # Output point with the specified type
        if out_type == "polar":
            p = self.shape.cart2polar(p[0], p[1], p[2])
        return p
        
    
    
class AreaRejE(EllipsoidBase):
    
    def __init__(self, shape, rng_type="MT", seed=12345):
        """! Initialiser for the area rejection (polar coordinates) ellipsoid sampler.
        @param shape : EllipsoidShape \n
            The ellipsoid whose surface is to be sampled.
        @param rng_type : str (optional) \n
            Specifies the random number generator to use. There are two options: \n
            - "MT"  : Mersenne twister (default setting).
            - "PCG" : PCG64.
        @param seed : int (optional) \n
            The seed of the random number generator, defaults to 12345.
        """
        self.initialise(shape, rng_type, seed)
        # Pre-calculate frequently used values
        if shape.is_spheroid() == True:
            a2 = shape.a_axis*shape.a_axis
            c2 = shape.c_axis*shape.c_axis
            self.a4 = a2*a2
            self.ac = a2*c2
            if self.shape.a_axis <= self.shape.c_axis: # i.e. if prolate or sphere
                self.M = 1.0 / (self.a_axis*self.c_axis)
            else: # i.e. if oblate
                sin_2_th = -a2 / (2.0*(c2-a2))
                sin_th_max = sqrt(sin_2_th)
                self.M = sin_th_max * sqrt(self.ac*sin_2_th + self.a4*(1.0 - sin_2_th))
                self.M = 1.0 / self.M
        else: # Triaxial
            a2 = shape.a_axis*shape.a_axis
            b2 = shape.b_axis*shape.b_axis
            c2 = shape.c_axis*shape.c_axis
            self.bc = b2*c2
            self.ac = a2*c2
            self.ab = a2*b2
            sin_th_max = sqrt( -b2 / (2.0*(c2-b2)) )
            sin_2_th = sin_th_max * sin_th_max
            self.M = sin_th_max * sqrt(self.ac*sin_2_th + self.ab*(1.0 - sin_2_th))
            self.M = 1.0 / self.M
        
    def random_surface_point(self, out_type="carts"):
        """! Generates a uniformly random point on the surface of an ellipsoid, 
        using the gradient rejection method.
        @param out_type : str (optional) \n
            Specifies the output type of the point.
            - "carts"  : Cartesian coordinates, \f$(x,y,z)\f$ (default setting).
            - "polar"  : Spherical polar coorinates, \f$(\theta, \phi)\f$.
        @return NumPy array \n
            Size of the output array depends on the input.
            - Cartesian coordinates of the generated surface point (if out_type = "carts").
            - Polar coordinates of the generated surface point (if out_type = "polar").
        """        
        accept = False
        if self.shape.is_spheroid() == True:
            ph = 2.0*pi*self.rng.random()
            sin_ph = sin(ph)
            sin_2_ph = sin_ph * sin_ph
            while accept == False:
                th = pi*self.rng.random()
                sin_th = sin(th)
                sin_2_th = sin_th*sin_th
                # Acceptance probability
                m = self.M * sin_th * sqrt(self.ac*sin_2_th + self.a4*(1.0 - sin_2_th))
                u = self.rng.random() # 'Roll dice'
                if u <= m:
                    accept = True
                self.no_attempts += 1
        else:
            while accept == False:
                th = pi*self.rng.random()
                ph = 2.0*pi*self.rng.random()
                sin_th = sin(th)
                sin_2_th = sin_th*sin_th
                sin_ph = sin(ph)
                sin_2_ph = sin_ph*sin_ph
                # Acceptance probability
                m = self.M * sin_th * sqrt( sin_2_th*(self.bc*(1.0-sin_2_ph) + self.ac*sin_2_ph) + self.ab*(1.0 - sin_2_th) )
                u = self.rng.random() # 'Roll dice'
                if u <= m:
                    accept = True
                self.no_attempts += 1
        
        # Output point with the specified type
        if out_type == "carts":
            p = self.shape.polar2cart(th, ph)
        else:
            p = array([0.0]*2)
            p[0] = th
            p[1] = ph

        return p
    
 
    
class RayIntersect(EllipsoidBase):
    
    def __init__(self, shape, rng_type="MT", seed=12345):
        """! Initialiser for the generic surface ellipsoid sampler.
        @param shape : EllipsoidShape \n
            The ellipsoid whose surface is to be sampled.
        @param rng_type : str (optional) \n
            Specifies the random number generator to use. There are two options: \n
            - "MT"  : Mersenne twister (default setting).
            - "PCG" : PCG64.
        @param seed : int (optional) \n
            The seed of the random number generator, defaults to 12345.
        """
        self.initialise(shape, rng_type, seed)
        self.sphere_sampler = Trig(1.0)
        self.sphere_sampler.set_rng(self.rng) # Make both sphere and ellipsoid samplers use the same rng object.
        # Pre-calculate frequently used values
        self.a2 = 1.0 / (shape.a_axis*shape.a_axis)
        self.b2 = 1.0 / (shape.b_axis*shape.b_axis)
        self.c2 = 1.0 / (shape.c_axis*shape.c_axis)
        self.R = shape.a_axis
        if shape.b_axis > self.R:
            self.R = shape.b_axis
        if shape.c_axis > self.R:
            self.R = shape.c_axis
        
    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of an ellipsoid, 
        using the generic surface sampler method of Detwiler et al.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """        
        accept = False
        while accept == False:
            # Generate sphere point, i.e. normal vector at point on a unit sphere
            n = self.sphere_sampler.random_surface_point() 
            # Find two orthogonal vectors in the tangent plane
            u = array([0.0, 0.0, 0.0])
            v = array([0.0, 0.0, 0.0])
            u[0] = 1.0
            u[1] = 1.0
            u[2] = -(n[0] + n[1])/n[2]
            v[0] = n[2] - n[1]*u[2]
            v[1] = n[0]*u[2] - n[2]
            v[2] = n[1] - n[0]
            u /= norm(u)
            v /= norm(v)
            # Generate point on disk tangent to bounding sphere at point r
            psi = 2.0*pi*self.rng.random()
            B = self.R*sqrt(self.rng.random()) # B^2 from uniform distribution in (0, R)
            b0 = cos(psi)
            b1 = B*sqrt(1.0 - b0*b0)
            b0 *= B
            if psi > pi:
                b1 *= -1.0 # If psi > pi, then sin(psi) is negative
            # Convert point in plane to 3D Cartesians
            p = self.R*n + b0*u + b1*v
            # Find intersections of line p0 + tp with the ellipsoid
            alpha = self.a2*n[0]*n[0] + self.b2*n[1]*n[1] + self.c2*n[2]*n[2]
            beta = -2.0 * (self.a2*n[0]*p[0] + self.b2*n[1]*p[1] + self.c2*n[2]*p[2])
            gamma = self.a2*p[0]*p[0] + self.b2*p[1]*p[1] + self.c2*p[2]*p[2] - 1.0
            t = beta*beta - 4.0*alpha*gamma
            if t < 0.0: # No intersections => try again
                accept = False
            else:
                j = 1.0 + round(self.rng.random()) # Random integer (1 or 2)
                if t > 1.0e-12: # Two intersections => select one of them at random
                    if j == 1:
                        t = 0.5 * (-beta + sqrt(t)) / alpha
                    else:
                        t = 0.5 * (-beta - sqrt(t)) / alpha
                elif fabs(t) < 1.0e-12 and j == 1: # One intersection and accepted
                    t = -0.5*beta/alpha
                accept = True
            self.no_attempts += 1
        return p - t*n
    
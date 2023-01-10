# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 11:38:46 2023

@author: Callum Marples

Random sampling of the surface of a sphere/ellipsoid.
"""

from math import sqrt, pi, cos, sin
from numpy import array
from numpy.linalg import norm
from numpy.random import Generator, PCG64, MT19937

### Sphere Samplers
class SphereBase:
    
    def __init__(self, radius=1.0, rng_type="MT", seed=12345):
        """! The initialiser for all sphere sampling classes.
        @param radius : float (optional) \n
            Radius of the sphere, defaults to 1.0.
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
        self.radius = radius
        self.no_attempts = 0
        
    
class SphereRejection(SphereBase):
    
    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of a sphere, 
        using the cubic rejection method.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """
        p = array([0.0, 0.0, 0.0])
        r_sq = 2.0
        # Generate a random point in a cube of length 2, accepting it if it 
        # lies within a unit sphere.
        while r_sq > 1.0:
            p[0] = 1.0 - 2.0*self.rng.random()
            p[1] = 1.0 - 2.0*self.rng.random()
            p[2] = 1.0 - 2.0*self.rng.random()
            r_sq = p[0]*p[0] + p[1]*p[1] + p[2]*p[2]
            self.no_attempts += 1
        # Scale the point so that it lies on the surface of the desired sphere.
        r_sq = 1.0 / sqrt(r_sq)
        p *= r_sq * self.radius
        return p
    

class SphereMarsaglia(SphereBase):
    
    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of a sphere, 
        using Marsaglia's method.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """
        p = array([0.0, 0.0, 0.0])
        r_sq = 2.0
        # Generate a random point in a cube of length 2, accepting it if it 
        # lies within a unit sphere.
        while r_sq >= 1.0:
            r1 = 1.0 - 2.0*self.rng.random()
            r2 = 1.0 - 2.0*self.rng.random()
            r_sq = r1*r1 + r2*r2
            self.no_attempts += 1
        p[2] = self.radius * (1.0 - 2.0*r_sq)
        r_sq = 2.0*sqrt(1.0 - r_sq)
        p[0] = self.radius*r1*r_sq;
        p[1] = self.radius*r2*r_sq;
        return p
    
    def set_rng(self, rng_in):
        self.rng = rng_in
    
class SphereCook(SphereBase):

    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of a sphere, 
        using Cook's method.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """
        p = array([0.0, 0.0, 0.0])
        s = 2.0
        while s >= 1.0:
            v1 = -1.0 + 2.0*self.rng.random()
            v2 = -1.0 + 2.0*self.rng.random()
            v3 = -1.0 + 2.0*self.rng.random()
            v4 = -1.0 + 2.0*self.rng.random()
            s = v1*v1 + v2*v2 + v3*v3 + v4*v4
            self.no_attempts += 1
        s = 1.0 / s
        p[0] = self.radius * 2.0 * s * (v2*v4 + v1*v3);
        p[1] = self.radius * 2.0 * s * (v3*v4 - v1*v2);
        p[2] = self.radius * s * (v1*v1 + v4*v4 - v2*v2 - v3*v3);
        return p;
            
class SphereGaussian(SphereBase):

    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of a sphere, 
        using Gaussian random numbers (with mean 0 and standard deviation 1).
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """            
        p = array([0.0, 0.0, 0.0])
        p[0] = self.rng.normal(0.0, 1.0)
        p[1] = self.rng.normal(0.0, 1.0)
        p[2] = self.rng.normal(0.0, 1.0)
        p /= norm(p)
        self.no_attempts += 1;
        return (self.radius * p);    
            
class SphereTrig(SphereBase):

    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of a sphere, 
        using the trigonometric method.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """            
        p = array([0.0, 0.0, 0.0])    
        u = 1.0 - 2.0*self.rng.random()
        ph = 2.0*pi*self.rng.random()
        p[2] = self.radius * u
        u = self.radius * sqrt(1.0 - u*u)
        cos_phi = cos(ph)
        p[0] = u * cos_phi
        p[1] = u * sqrt(1.0 - cos_phi*cos_phi)
        if ph > pi:
            p[1] *= -1.0
        self.no_attempts += 1
        return p
    
    def set_rng(self, rng_in):
        self.rng = rng_in
    
class SphereAreaReject(SphereBase):

    def random_surface_point(self):
        """! Generates a uniformly random point on the surface of a sphere, 
        using the area element rejection method.
        @return NumPy array \n
            Cartesian coordinates of the generated surface point.
        """        
        p = array([0.0, 0.0, 0.0]) 
        accept = False
        ph = 2.0*pi*self.rng.random()
        while accept == False:
            th = pi*self.rng.random()
            sin_th = sin(th)      # Acceptance probability
            u = self.rng.random() # 'Roll dice'
            if u <= sin_th:
                accept = True
            self.no_attempts += 1
        p[2] = self.radius*sqrt(1.0 - sin_th*sin_th)
        if th > 0.5*pi:
            p[2] *= -1.0
        sin_th *= self.radius
        cos_ph = cos(ph)
        p[0] = sin_th * cos_ph
        p[1] = sin_th * sqrt(1.0 - cos_ph*cos_ph)
        if ph > pi:
            p[1] *= -1.0
        return p
        
### Ellipsoid Samplers

class EllipsoidBase:
    
    #def __init__(self, a, b, c, rng_type="MT", seed=12345):
    #    self.initialise(a, b, c, rng_type, seed)
    
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


class EllipsoidNaive(EllipsoidBase):
    
    def __init__(self, shape, rng_type="MT", seed=12345):
        """! Initialiser for the naive sampler.
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
        self.initialise(shape, rng_type, seed)
        self.sphere_sampler = SphereMarsaglia(1.0)
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
        
        
        









# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:32:28 2023

@author: Callum Marples

Examples of using the sampling routines.
"""

import leod
import os

os.chdir("..")

### Sphere Samplers

## Cubic rejection
S1 = leod.sampling.SphereRejection()  # Create sampler
p1 = S1.random_surface_point()        # Use sampler to generate random point

## Marsaglia's method
S2 = leod.sampling.SphereMarsaglia()  # Create sampler
p2 = S2.random_surface_point()        # Use sampler to generate random point

## Cook's method
S3 = leod.sampling.SphereCook()       # Create sampler
p3 = S3.random_surface_point()        # Use sampler to generate random point

## Gaussian method
S4 = leod.sampling.SphereGaussian()   # Create sampler
p4 = S4.random_surface_point()        # Use sampler to generate random point

## Trig method
S5 = leod.sampling.SphereTrig()       # Create sampler
p5 = S5.random_surface_point()        # Use sampler to generate random point

## Area rejection method
S6 = leod.sampling.SphereAreaReject() # Create sampler
p6 = S6.random_surface_point()        # Use sampler to generate random point



### Ellipsoid Samplers

# Define shape
shape = leod.shape.EllipsoidShape(3.0, 2.0, 1.0)

## Naive scaling
E1 = leod.sampling.EllipsoidNaive(shape)
r1 = E1.random_surface_point()


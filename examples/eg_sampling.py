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
S1 = leod.sampling.sphere.CubicReject()      # Create sampler
p1 = S1.random_surface_point()        # Use sampler to generate random point

## Marsaglia's method
S2 = leod.sampling.sphere.Marsaglia()        
p2 = S2.random_surface_point()       

## Cook's method
S3 = leod.sampling.sphere.Cook()       
p3 = S3.random_surface_point()        

## Gaussian method
S4 = leod.sampling.sphere.Gaussian()   
p4 = S4.random_surface_point()        

## Trig method
S5 = leod.sampling.sphere.Trig()       
p5 = S5.random_surface_point()        

## Area rejection method
S6 = leod.sampling.sphere.AreaReject() 
p6 = S6.random_surface_point()        



### Ellipsoid Samplers

# Define shape
shape = leod.shape.EllipsoidShape(3.0, 2.0, 1.0)
shape2 = leod.shape.EllipsoidShape(2.0, 2.0, 1.0)

## Naive scaling
E1 = leod.sampling.ellipsoid.NaiveScaling(shape)
r1 = E1.random_surface_point()

## Gradient rejection (with Marsaglia's method)
E2 = leod.sampling.ellipsoid.GradRej(shape)
r2 = E2.random_surface_point()

## Gradient rejection (with the 'Trig method')
E3 = leod.sampling.ellipsoid.GradRej(shape, sphere_type="Trig")
r3 = E3.random_surface_point()

## Area rejection (polar)
E4 = leod.sampling.ellipsoid.AreaRejPol(shape)
r4 = E4.random_surface_point()
r4p = E4.random_surface_point("polar")

## Area rejection (polar) on a spheroid
E5 = leod.sampling.ellipsoid.AreaRejPol(shape2)
r5 = E5.random_surface_point()
r5p = E5.random_surface_point("polar")

## Area rejection (Mercator)
E6 = leod.sampling.ellipsoid.AreaRejMer(shape)
r6 = E6.random_surface_point()

## Area rejection (Mercator) on a spheroid
E7 = leod.sampling.ellipsoid.AreaRejMer(shape2)
r7 = E7.random_surface_point()

## Generic sampler
E8 = leod.sampling.ellipsoid.Generic(shape)
r8 = E8.random_surface_point()


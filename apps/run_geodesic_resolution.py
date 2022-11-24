# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:18:06 2022

@author: Callum Marples

Calculate distance for a particular example, at a range of resolutions.

Example: (theta, phi) = (50, 60) -> (90, 0)
(above angles in degrees)

Five different shapes are used and for each, an appropriate alternate method
is used for comparison.
"""

import leod
import math
import numpy as np
import csv
import sys
import os

modulename = 'leod.spheroid_geodesics'
if modulename in sys.modules:
    geo_flag = True
else:
    geo_flag = False

os.chdir("..")

# Define shapes.
shape = [ leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 1.0),   # Sphere.
          leod.ellipsoid_shape.EllipsoidShape(6378137.0, 6378137.0, 6356752.3142),  # Earth (WGS84 ellipsoid).
          leod.ellipsoid_shape.EllipsoidShape(2.0, 2.0, 1.0),   # Oblate.
          leod.ellipsoid_shape.EllipsoidShape(1.0, 1.0, 2.0),   # Prolate.
          leod.ellipsoid_shape.EllipsoidShape(3.0, 2.0, 1.0) ]  # Triaxial.

# Normalise shapes so that (abc) = 1.
# Unit sphere already normalised, do not want to do this for the WGS84 example.
shape[2].normalise()
shape[3].normalise()
shape[4].normalise()

# Define start and end points of the geodesic example.
start_point = [50.0, 60.0]
end_point = [90.0, 0.0]
conv = math.pi / 180.0 # Degrees to radians.

# Number of vertices (logarithmic scale)
log_n = np.linspace(2, 6, 25)
no_vertices_float = np.power(10.0, log_n)
no_vertices = [0] * len(no_vertices_float)
for i in range(len(no_vertices)):
    no_vertices[i] = int(round(no_vertices_float[i]))

# Find mesh parameters from the number of vertices.
no_vertices_polar = [0] * len(no_vertices)
no_vertices_icosahedral = [0] * len(no_vertices)
no_divisions = [0] * len(no_vertices)
no_phi = [0] * len(no_vertices)
no_theta = [0] * len(no_vertices)
for i in range(len(no_vertices)):
    # Solve quadratic 10*n_d^2 + 2 - n_v = 0 for n_d, the number of face
    # divisions in the icosahedral triangulation.
    no_divisions[i] = round(math.sqrt(0.1*(no_vertices[i] - 2)))
    no_vertices_icosahedral[i] = 10*no_divisions[i]**2 + 2
    # Solve quadratic 0.5*n_phi^2 - n_phi + 2 - n_v = 0 for n_phi, the number
    # of phi values in the theta-phi grid.
    no_phi[i] = round(1 + math.sqrt(2*no_vertices[i] - 3))
    if no_phi[i] % 2 == 1:
        # Make no_phi even by adding one.
        no_phi[i] += 1
    no_theta[i] = round(0.5*no_phi[i]) + 1
    no_vertices_polar[i] = no_phi[i]*(no_theta[i]-2) + 2


            
        
        
        
        
        
# Initialise result lists
    
    

# Compute the test example distance on each shape in turn, for each resolution
# and mesh type.
#for i in range(len(shape)):

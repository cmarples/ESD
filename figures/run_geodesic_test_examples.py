# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 14:58:35 2022

@author: Callum Marples

Run a set of test examples using the fast marching method, along with any other
methods appropriate for the given shape.
"""

# Imports

from leod.ellipsoid_shape import EllipsoidShape
from leod.fmm_polar_graph import generate_polar_graph
from leod.fmm_precalculation import precalculate_grid
from leod.fmm_callers import calculate_pair_distance

# Test examples
# Element 0 of each sub-list is the theta coordinate, element 1 is the phi coordinate.
# All angular coordinates are in degrees.
start_points = [ [0.0, 0.0],
                 [0.0, 0.0],
                 [0.0, 0.0],
                 [90.0, 0.0],
                 [90.0, 0.0],
                 [50.0, 60.0],
                 [10.0, 25.0],
                 [10.0, 25.0],
                 [70.0, 80.0],
                 [120.0, 40.0]] 

end_points = [ [180.0, 0.0],
               [90.0, 0.0],
               [90.0, 90.0],
               [50.0, 60.0],
               [130.0, 60.0],
               [90.0, 0.0],
               [20.0, 35.0],
               [10.0, 190.0],
               [40.0, 340.0],
               [60.0, 220.0] ] 

shape = [ EllipsoidShape(1.0, 1.0, 1.0),   # Sphere
          EllipsoidShape(6378137.0, 6378137.0, 6356752.3142),  # Earth (WGS84 ellipsoid)
          EllipsoidShape(2.0, 2.0, 1.0),   # Oblate
          EllipsoidShape(1.0, 1.0, 2.0),   # Prolate
          EllipsoidShape(3.0, 2.0, 1.0) ]  # Triaxial
    
d8_fmm = []
obtuse8 = []

# Define theta-phi grid parameters
no_theta = 181
no_phi = 360

# Compute test examples on each shape in turn
#for i in range(len(shape)):
for i in range(5):
    
    # Initialise lists
    d8_fmm.append([])
    obtuse8.append([])
    
    # Scale ellipsoid axes in each shape so that abc = 1
    # i.e. so each have the volume of a unit sphere (4*pi/3)
    if i != 2: # Do not normalise the Earth example
        shape[i].normalise()
        
    # Generate theta-phi (polar) grid and precalculate distances and face angles
    vertex, polar_grid = generate_polar_graph(shape[i], no_theta, no_phi, is_connect_8=True, is_Dijkstra=False)
    max_angle = precalculate_grid(vertex)

    # Call the pair distance fmm routine for each start-end pair
    for j in range(len(start_points)):
        d, fmm = calculate_pair_distance(shape[i], vertex, start_points[j], end_points[j],
                                        1, graph_type="polar", grid=polar_grid, is_radians=False)
        d8_fmm[i].append(d)
        obtuse8[i].append(fmm.no_obtuse_path)
    
    
    
    
       



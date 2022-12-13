# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 16:51:22 2022

@author: Cal
"""

import leod
import math
import copy
import os

os.chdir("..")

dijkstra_flag = True

conv = math.pi / 180.0 # Degrees to radians.
start_point = [50.0*conv, 60.0*conv]
end_point = [90.0*conv, 0.0*conv]

### Sphere

# Taxicab
sphere = leod.shape.EllipsoidShape(1.0)
d_taxi_sphere = leod.geo.taxicab.sphere_td(sphere.a_axis, start_point, end_point)
print("")
print("Sphere:             taxicab distance =", d_taxi_sphere)

### Spheroid
sph = leod.shape.EllipsoidShape(2.0, 2.0, 1.0)
sph.normalise()
d_taxi_spheroid = leod.geo.taxicab.spheroid_td(sph.a_axis, sph.c_axis, start_point, end_point)
print("Oblate spheroid:    taxicab distance =", d_taxi_spheroid)

### Triaxial
tri = leod.shape.EllipsoidShape(3.0, 2.0, 1.0)
tri.normalise()
d_taxi_triaxial = leod.geo.taxicab.triaxial_td(tri.a_axis, tri.b_axis, tri.c_axis, start_point, end_point)
print("Triaxial ellipsoid: taxicab distance =", d_taxi_triaxial)
print("")

if dijkstra_flag == True:
    
    no_theta = 91
    no_phi = 180
    
    grid_sphere = leod.fmm.grid_pol.gen_pol_grid(no_theta, no_phi, is_connect_8=False, is_generic=True)
    grid_sph = copy.deepcopy(grid_sphere)
    grid_tri = copy.deepcopy(grid_sphere)
    
    leod.fmm.grid_pol.pre_pol_grid(grid_sphere, sphere)
    leod.fmm.grid_pol.pre_pol_grid(grid_sph, sph)
    leod.fmm.grid_pol.pre_pol_grid(grid_tri, tri)
    
    d_sphere, fmm1 = leod.fmm.callers.distance_pair(sphere, grid_sphere, start_point, end_point,
                                                   is_radians=True, is_dijkstra=True)
    d_sph, fmm2 = leod.fmm.callers.distance_pair(sph, grid_sph, start_point, end_point,
                                                   is_radians=True, is_dijkstra=True)
    d_tri, fmm3 = leod.fmm.callers.distance_pair(tri, grid_tri, start_point, end_point,
                                                   is_radians=True, is_dijkstra=True)
    print("Sphere:             Dijkstra-4 distance =", d_sphere)
    print("Oblate spheroid:    Dijkstra-4 distance =", d_sph)
    print("Triaxial ellipsoid: Dijkstra-4 distance =", d_tri)
    print("")
    
    
    
    
    
    
    
    
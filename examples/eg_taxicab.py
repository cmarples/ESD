"""
@brief Example script for sphere, spheroid and ellipsoid taxicab distances.
@file eg_taxicab.py
@author Callum Marples

- Created on 12/12/2022. 
- Last modified on 17/01/2023.
"""

import leod
import math
import os

os.chdir("..")

dijkstra_flag = True

start_point = [50.0, 60.0]
end_point = [90.0, 0.0]

### Taxicab distances
# Sphere
sphere = leod.shape.EllipsoidShape(1.0, 1.0, 1.0)
d_taxi_sphere = leod.geo.taxicab.sphere_tcd(sphere.a_axis, start_point, end_point, out_flag=True, is_radians=False)
print("")
print("Sphere:             taxicab distance =", d_taxi_sphere[0])

# Spheroid
sph = leod.shape.EllipsoidShape(2.0, 2.0, 1.0)
sph.normalise()
d_taxi_spheroid = leod.geo.taxicab.spheroid_tcd(sph.a_axis, sph.c_axis, start_point, end_point, out_flag=True, is_radians=False)
print("Oblate spheroid:    taxicab distance =", d_taxi_spheroid[0])

# Triaxial
tri = leod.shape.EllipsoidShape(3.0, 2.0, 1.0)
tri.normalise()
d_taxi_tri = leod.geo.taxicab.triaxial_tcd(tri, start_point, end_point, out_flag=True, is_radians=False)
print("Triaxial ellipsoid: taxicab distance =", d_taxi_tri[0])
print("")

### Dijkstra's algorithm (optional: for comparison)
if dijkstra_flag == True:
    
    no_theta = 91
    no_phi = 180
    
    grid_sphere = leod.fmm.mesh_pol.gen_pol_mesh(no_theta, no_phi, shape=sphere, is_connect_8=False)
    grid_sph = leod.fmm.mesh_pol.gen_pol_mesh(no_theta, no_phi, shape=sph, is_connect_8=False)
    grid_tri = leod.fmm.mesh_pol.gen_pol_mesh(no_theta, no_phi, shape=tri, is_connect_8=False)
    
    d_sphere, fmm1 = leod.fmm.callers.distance_pair(sphere, grid_sphere, start_point, end_point,
                                                   is_radians=False, is_dijkstra=True)
    d_sph, fmm2 = leod.fmm.callers.distance_pair(sph, grid_sph, start_point, end_point,
                                                   is_radians=False, is_dijkstra=True)
    d_tri, fmm3 = leod.fmm.callers.distance_pair(tri, grid_tri, start_point, end_point,
                                                   is_radians=False, is_dijkstra=True)
    print("Sphere:             Dijkstra-4 distance =", d_sphere)
    print("Oblate spheroid:    Dijkstra-4 distance =", d_sph)
    print("Triaxial ellipsoid: Dijkstra-4 distance =", d_tri)
    print("")
    
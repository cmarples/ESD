"""
@brief Example script for mesh generation.
@file eg_mesh_generation.py
@author Callum Marples

- Created on 26/01/2023. 
- Last modified on 26/01/2023.
"""

from esd.shape import EllipsoidShape
from esd.fmm.mesh_pol import gen_pol_mesh
from esd.fmm.mesh_ico import gen_ico_mesh
from math import sqrt, acos, pi
import time

### Grid size
# Set the parameters of each grid so their resolutions are as close as possible. 
no_vertices = 4000
# Solve quadratic 10*n_d^2 + 2 - n_v = 0 for n_d, the number of face
# divisions in the icosahedral triangulation.
no_divisions = round(sqrt(0.1*(no_vertices - 2)))
no_vertices_icosahedral = 10*no_divisions**2 + 2
# Solve quadratic 0.5*n_phi^2 - n_phi + 2 - n_v = 0 for n_phi, the number
# of phi values in the theta-phi grid.
no_phi = round(1 + sqrt(2*no_vertices - 3))
if no_phi % 2 == 1:
    # Make no_phi even by adding one.
    no_phi += 1
no_theta = round(0.5*no_phi) + 1
no_vertices_polar = no_phi*(no_theta-2) + 2

### Define ellipsoid
shape = EllipsoidShape(3.0, 2.0, 1.0)
shape.normalise()

#### Polar grid with 4 neighbours
print("")
print("Run the fast marching method on a 4-neighbour polar mesh.")
tic = time.perf_counter()

mesh1 = gen_pol_mesh(no_theta, no_phi, shape, is_connect_8=False)

toc = time.perf_counter()
t1 = toc - tic
print("Generation run-time =", t1)
print("Number of vertices  =", mesh1.no_vertices)
print("Minimum edge length =", mesh1.min_edge)
print("Maximum edge length =", mesh1.max_edge)
print("No. obtuse angles   =", mesh1.no_obtuse)
print("Minimum face angle  =", acos(mesh1.min_angle)*180.0/pi, "degrees") # Output as angles in degrees
print("Maximum face angle  =", acos(mesh1.max_angle)*180.0/pi, "degrees")

### Polar grid with 8 neighbours
print("")
print("Run the fast marching method on an 8-neighbour polar mesh.")
tic = time.perf_counter()

mesh2 = gen_pol_mesh(no_theta, no_phi, shape, is_connect_8=True)

toc = time.perf_counter()
t2 = toc - tic
print("Generation run-time =", t2)
print("Number of vertices  =", mesh2.no_vertices)
print("Minimum edge length =", mesh2.min_edge)
print("Maximum edge length =", mesh2.max_edge)
print("No. obtuse angles   =", mesh2.no_obtuse)
print("Minimum face angle  =", acos(mesh2.min_angle)*180.0/pi, "degrees") # Output as angles in degrees
print("Maximum face angle  =", acos(mesh2.max_angle)*180.0/pi, "degrees")

### Icosahedral grid
print("")
print("Run the fast marching method on an icosahedral mesh.")
tic = time.perf_counter()

mesh3 = gen_ico_mesh(no_divisions, shape, is_split=False)

toc = time.perf_counter()
t3 = toc - tic
print("Generation run-time =", t3)
print("Number of vertices  =", mesh3.no_vertices)
print("Minimum edge length =", mesh3.min_edge)
print("Maximum edge length =", mesh3.max_edge)
print("No. obtuse angles   =", mesh3.no_obtuse)
print("Minimum face angle  =", acos(mesh3.min_angle)*180.0/pi, "degrees") # Output as angles in degrees
print("Maximum face angle  =", acos(mesh3.max_angle)*180.0/pi, "degrees")

### Icosahedral grid with splitting
print("")
print("Run the fast marching method on a split icosahedral mesh.")
tic = time.perf_counter()

mesh4 = gen_ico_mesh(no_divisions, shape, is_split=True)

toc = time.perf_counter()
t4 = toc - tic
print("Generation run-time =", t4)
print("Number of vertices  =", mesh4.no_vertices)
print("Minimum edge length =", mesh4.min_edge)
print("Maximum edge length =", mesh4.max_edge)
print("No. obtuse angles   =", mesh4.no_obtuse)
print("Minimum face angle  =", acos(mesh4.min_angle)*180.0/pi, "degrees") # Output as angles in degrees
print("Maximum face angle  =", acos(mesh4.max_angle)*180.0/pi, "degrees")

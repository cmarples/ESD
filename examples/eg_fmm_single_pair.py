"""
@brief Example script for mesh generation and the single pair fast marching distance routine.
@file eg_fmm_single_pair.py
@author Callum Marples

To use Dijkstra's algorithm, set is_dijkstra=True in the distance_pair function.

- Created on 16/01/2023. 
- Last modified on 17/01/2023.
"""

from leod.shape import EllipsoidShape
from leod.fmm.mesh_pol import gen_pol_mesh
from leod.fmm.mesh_ico import gen_ico_mesh
from leod.fmm.callers import distance_pair
import math
import time

### Grid size
# Set the parameters of each grid so their resolutions are as close as possible. 
no_vertices = 4000
# Solve quadratic 10*n_d^2 + 2 - n_v = 0 for n_d, the number of face
# divisions in the icosahedral triangulation.
no_divisions = round(math.sqrt(0.1*(no_vertices - 2)))
no_vertices_icosahedral = 10*no_divisions**2 + 2
# Solve quadratic 0.5*n_phi^2 - n_phi + 2 - n_v = 0 for n_phi, the number
# of phi values in the theta-phi grid.
no_phi = round(1 + math.sqrt(2*no_vertices - 3))
if no_phi % 2 == 1:
    # Make no_phi even by adding one.
    no_phi += 1
no_theta = round(0.5*no_phi) + 1
no_vertices_polar = no_phi*(no_theta-2) + 2

### Start and end points
start_point = [90.0, 0.0]
end_point = [50.0, 60.0]

### Define ellipsoid
shape = EllipsoidShape(1.0, 1.0, 1.0)
shape.normalise()

#### Polar grid with 4 neighbours
print("")
print("Run the fast marching method on a 4-neighbour polar mesh.")
tic = time.perf_counter()

mesh1 = gen_pol_mesh(no_theta, no_phi, shape, is_connect_8=False)
d1, fmm1 = distance_pair(shape, mesh1, start_point, end_point, is_dijkstra=False)

toc = time.perf_counter()
t1 = toc - tic
print("Distance =", d1)
print("Run-time =", t1)

### Polar grid with 8 neighbours
print("")
print("Run the fast marching method on an 8-neighbour polar mesh.")
tic = time.perf_counter()

mesh2 = gen_pol_mesh(no_theta, no_phi, shape, is_connect_8=True)
d2, fmm2 = distance_pair(shape, mesh2, start_point, end_point, is_dijkstra=False)

toc = time.perf_counter()
t2 = toc - tic
print("Distance =", d2)
print("Run-time =", t2)

### Icosahedral grid
print("")
print("Run the fast marching method on an icosahedral mesh.")
tic = time.perf_counter()

mesh3 = gen_ico_mesh(no_divisions, shape, is_split=False)
d3, fmm3 = distance_pair(shape, mesh3, start_point, end_point, is_dijkstra=False)

toc = time.perf_counter()
t3 = toc - tic
print("Distance =", d3)
print("Run-time =", t3)

### Icosahedral grid with splitting
print("")
print("Run the fast marching method on a split icosahedral mesh.")
tic = time.perf_counter()

mesh4 = gen_ico_mesh(no_divisions, shape, is_split=True)
d4, fmm4 = distance_pair(shape, mesh4, start_point, end_point, is_dijkstra=False)

toc = time.perf_counter()
t4 = toc - tic
print("Distance =", d4)
print("Run-time =", t4)

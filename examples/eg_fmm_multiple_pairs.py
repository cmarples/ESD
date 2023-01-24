"""
@brief Example script for the multiple pair fast marching distance routine.
@file eg_fmm_multiple_pairs.py
@author Callum Marples

To use Dijkstra's algorithm, set is_dijkstra=True in the distance_pair function.

- Created on 17/01/2023. 
- Last modified on 17/01/2023.
"""

from esd.shape import EllipsoidShape
from esd.fmm.mesh_pol import gen_pol_mesh
from esd.fmm.callers import distance_multiple
import time

### Grid size
no_theta = 51
no_phi = 100

### Start point
start_point = [90.0, 0.0]

### End points
end_point = [[50.0, 60.0], [130.0, 60.0], [0.0, 0.0], [180.0, 0.0], [90.0, 90.0]]

### Define ellipsoid
shape = EllipsoidShape(3.0, 2.0, 1.0)
shape.normalise()

print("")
print("Run the fast marching method on an 8-neighbour polar mesh, with 5 end points.")
tic = time.perf_counter()

### Define a Polar mesh with 8 neighbours
mesh = gen_pol_mesh(no_theta, no_phi, shape, is_connect_8=True)

### Run the fast marching method
d, fmm = distance_multiple(shape, mesh, start_point, end_point, is_dijkstra=False)

toc = time.perf_counter()
t = toc - tic
print("Distance =", d)
print("Run-time =", t)


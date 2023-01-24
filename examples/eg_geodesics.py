"""
@brief Example script for sphere, spheroid and ellipsoid geodesics (non-wavefront).
@file eg_geodesics.py
@author Callum Marples

- Created on 17/01/2023. 
- Last modified on 17/01/2023.
"""

from esd.shape import EllipsoidShape
from esd.geo.sphere import gc_dist
from esd.geo.spheroid import geo_dist
from esd.geo.triaxial import bvm_dist
import os

os.chdir("..")

### Start and end points
start_point = [90.0, 0.0]
end_point = [50.0, 60.0]

### Sphere
r = 1.0 # Radius
# Use the great circle derived formula
d0 = gc_dist(r, start_point, end_point, is_radians=False)
print("Sphere : d = ", d0)

### Spheroid
sph = EllipsoidShape(2.0, 2.0, 1.0)
sph.normalise()
# Use GeographicLib
d1 = geo_dist(sph, start_point, end_point, is_radians=False)
print("Spheroid (Geo): d = ", d1)
# Use the boundary value method
d2 = bvm_dist(sph, start_point, end_point, is_radians=False, tol=1e-12, Jacobi=False, n=5000)
print("Spheroid (BVM): d = ", d2[0])

### Triaxial
tri = EllipsoidShape(3.0, 2.0, 1.0)
tri.normalise()
# Use the boundary value method
d3 = bvm_dist(tri, start_point, end_point, is_radians=False, tol=1e-12, Jacobi=False, n=5000)
print("Triaixal (BVM): d = ", d3[0])






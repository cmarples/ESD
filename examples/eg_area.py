"""
@brief Example script for surface areas of a sphere/spheroid/ellipsoid.
@file eg_area.py
@author Callum Marples

- Created on 12/01/2023. 
- Last modified on 17/01/2023.
"""

from math import pi
import leod
import os

os.chdir("..")

# Define shapes
sphere = leod.shape.EllipsoidShape(1.0, 1.0, 1.0)
sph = leod.shape.EllipsoidShape(2.0, 2.0, 1.0)
ell = leod.shape.EllipsoidShape(3.0, 2.0, 1.0)

### Full Surface Area

# Sphere
a1 = leod.area.full.sphere(1.0)

# Spheroid
a2 = leod.area.full.spheroid(sph)

# Thomsen
a3 = leod.area.full.thomsen(ell)

# Legendre
a4 = leod.area.full.legendre(ell)

# Numerical (Romberg)
a5 = leod.area.full.numerical(ell)

### Patch Areas

# Sphere patch
b1 = leod.area.patch.sphere(1.0, 0.5*pi, 0.6*pi, 0.2*pi, 0.3*pi)

# Spheroid band
b2 = leod.area.patch.spheroid_band(sph, 0.5*pi, 0.6*pi)

# Spheroid patch
b3 = leod.area.patch.spheroid(sph, 0.5*pi, 0.6*pi, 0.2*pi, 0.3*pi)

# Ellipsoid patch
b4 = leod.area.patch.ellipsoid(ell, 0.5*pi, 0.6*pi, 0.2*pi, 0.3*pi)
b5 = leod.area.patch.ellipsoid(sphere, 0.5*pi, 0.6*pi, 0.2*pi, 0.3*pi)
b6 = leod.area.patch.ellipsoid(sph, 0.5*pi, 0.6*pi, 0.2*pi, 0.3*pi)


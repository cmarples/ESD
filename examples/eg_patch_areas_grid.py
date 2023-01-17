"""
@brief Example script for surface areas of a sphere/spheroid/ellipsoid.
@file eg_patch_areas_grid.py
@author Callum Marples

- Created on 13/01/2023. 
- Last modified on 17/01/2023.
"""

import leod

from math import pi
import numpy as np
import os

os.chdir("..")

### Define shape and grid
ell = leod.shape.EllipsoidShape(3.0, 2.0, 1.0)
grid = leod.grid.Grid(181, 360)

### Compute patch areas using the patch.grid routine
p = leod.area.patch.grid(ell, grid)

# Does the combined areas of the patches agree with routines for full area?
a0 = np.sum(p) 
# Thomsen
a1 = leod.area.full.thomsen(ell)
# Legendre
a2 = leod.area.full.legendre(ell)
# Numerical (Romberg)
a3 = leod.area.full.numerical(ell)

# North pole patch area (should match p[0])
an = 4.0 * leod.area.patch.ellipsoid(ell, 0.0, 0.5*grid.delta_theta, 0.0, 0.5*pi)
# Pole adjacent patch area (should be half of p[1])
am = leod.area.patch.ellipsoid(ell, 0.5*grid.delta_theta, 1.5*grid.delta_theta, 0.0, 0.5*grid.delta_phi)
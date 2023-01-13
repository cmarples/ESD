# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 12:05:46 2023

@author: Callum Marples

Examples of using the area routines.
"""

import leod

from math import pi
import numpy as np
import os

os.chdir("..")

# Define shapes
ell = leod.shape.EllipsoidShape(3.0, 2.0, 1.0)

# Ellipsoid grid
grid = leod.grid.Grid(181, 360)

p = leod.area.patch.grid(ell, grid)
a0 = np.sum(p)

# Thomsen
a1 = leod.area.full.thomsen(ell)
# Legendre
a2 = leod.area.full.legendre(ell)
# Numerical (Romberg)
a3 = leod.area.full.numerical(ell)

# North pole patch area
an = 4.0 * leod.area.patch.ellipsoid(ell, 0.0, 0.5*grid.delta_theta, 0.0, 0.5*pi)

# Pole adjacent patch area
am = leod.area.patch.ellipsoid(ell, 0.5*grid.delta_theta, 1.5*grid.delta_theta, 0.0, 0.5*grid.delta_phi)
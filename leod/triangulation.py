# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 12:31:45 2022

@author: Callum Marples
"""

import math
import numpy as np

# Sphere triangulation using a geodesic polyhedron

# Icosahedron coordinates
tau = 0.5 * (1.0 + math.sqrt(5.0)) # Golden ratio
den = math.sqrt(tau*tau + 1.0)     # Divide by this to get unit circumradius

ico_vertex = [ [0, 1, tau], [0, -1, tau], [0, 1, -tau], [0, -1, -tau],
               [1, tau, 0], [-1, tau, 0], [1, -tau, 0], [-1, -tau, 0],
               [tau, 0, 1], [-tau, 0, 1], [tau, 0, -1], [-tau, 0, -1] ]
ico_vertex = np.array(ico_vertex) / den


# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:25:53 2022

@author: Callum Marples

Vertex class for use in the Fast Marching Method.
"""

import numpy as np

# This class contains information for a vertex in the grid used in the fast
# marching method.
class FmmVertex:
    def __init__(self, index, carts=np.zeros(3), is_endpoint=False):
        self.index = index
        self.is_endpoint = is_endpoint
        self.carts = carts
        self.neighbour = []
        self.neighbour_distance = []
        self.face = []
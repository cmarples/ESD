# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:25:53 2022

@author: Callum Marples

Vertex class for use in the Fast Marching Method.
"""

import math
import numpy as np

# This class contains information for a vertex in the grid used in the fast
# marching method.
class FmmVertex:
    def __init__(self, index, carts=np.zeros(3), is_endpoint=False):
        self.index = index
        self.is_endpoint = is_endpoint
        self.carts = carts
        self.neighbour = []
        self.face = []
        self.distance_to_vertex = {}



# Get the Euclidean distance between vertices i and j, updating the distance
# lists of both.
def get_distance(vertex, i, j):
    if j in vertex[i].distance_to_vertex:
        return vertex[i].distance_to_vertex[j]
    else:
        # Calculate distance if it has not already been found
        #vertex[i].distance_to_vertex[j] = math.sqrt( np.sum(vertex[i].carts - vertex[j].carts)**2 )
        vertex[i].distance_to_vertex[j] = math.sqrt( (vertex[i].carts[0]-vertex[j].carts[0])**2.0 +
                                                     (vertex[i].carts[1]-vertex[j].carts[1])**2.0 +
                                                     (vertex[i].carts[2]-vertex[j].carts[2])**2.0 )
        if i in vertex[j].neighbour:
            vertex[j].distance_to_vertex[i] = vertex[i].distance_to_vertex[j]
        return vertex[i].distance_to_vertex[j]        
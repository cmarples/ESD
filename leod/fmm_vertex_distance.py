# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:39:56 2022

@author: Callum Marples

Distance between two vertices
"""

import math
import numpy as np

# Get the Euclidean distance between vertices i and j, updating the distance
# lists of both.
def get_distance(vertex, i, j):
    if j in vertex[i].distance_to_vertex:
        return vertex[i].neighbour_distance[j]
    else:
        # Calculate distance if it has not already been found
        vertex[i].distance_to_vertex[j] = math.sqrt( np.sum(vertex[i].carts - vertex[j].carts)**2 )
        if i in j.neighbour:
            vertex[j].distance_to_vertex[i] = vertex[i].distance_to_vertex[j]
        return vertex[i].neighbour_distance[j]
        
    
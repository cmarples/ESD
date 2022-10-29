# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 12:31:45 2022

@author: Callum Marples
"""

import math
import numpy as np

# Check that each triangle is valid by calculating each vertex angle
def check_triangles(vertex):
    no_obtuse = 0
    # Check each vertex in turn
    for i in range(len(vertex)):
        for face_no in range(len(vertex[i].face)):
            # Triangle vertices
            j, k = vertex[i].face[face_no]
            # Edge vectors
            w1 = vertex[j].carts - vertex[i].carts
            w2 = vertex[k].carts - vertex[i].carts
            # Need angle between vectors w1 and w2 to be <= 90 degrees
            cos_alpha = np.dot(w1, w2)
            if cos_alpha < 0:
                no_obtuse += 1  
    return no_obtuse


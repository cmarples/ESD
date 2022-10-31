# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 12:31:45 2022

@author: Callum Marples
"""

import math
import numpy as np

# Find closest vertex to input point
def find_closest_vertex(vertex, p, v=0):

    d2 = {}
    is_min_found = False
    while is_min_found == False:
        if v not in d2.keys():
            vec = vertex[v].carts - p
            d2[v] = np.dot(vec, vec)
        u = v # Index of locally closest vertex
        d_min = d2[v]
        for j in vertex[v].distance_to_neighbour.keys():
            if j not in d2.keys():
                vec = vertex[j].carts - p
                d2[j] = np.dot(vec, vec)
            if d2[j] < d_min:
                d_min = d2[j]
                u = j
        if u == v:
            is_min_found = True
        else:
            v = u
    
    return v
    
    
    '''
    d_min = math.inf
    v_min = 0
    for i in range(len(vertex)):
        vec = vertex[i].carts - p
        d2 = np.dot(vec, vec)
        if d2 < d_min:
            d_min = d2
            v_min = i
    return v_min
    '''

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


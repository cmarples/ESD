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
        for j in vertex[v].neighbour.keys():
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
        vertex[i].face_dot = []
        for face_no in range(len(vertex[i].face)):
            # Triangle vertices
            j, k = vertex[i].face[face_no]
            # Edge vectors
            w1 = vertex[j].carts - vertex[i].carts
            w2 = vertex[k].carts - vertex[i].carts
            # Need angle between vectors w1 and w2 to be <= 90 degrees
            cos_alpha = np.dot(w1, w2)
            vertex[i].face_dot.append(cos_alpha)
            if cos_alpha < 0:
                no_obtuse += 1  
    return no_obtuse


# Find the face containing point p
def find_face(vertex, p, q):
    for i, j in vertex[q].face:
        u = vertex[i].carts - vertex[q].carts
        v = vertex[j].carts - vertex[q].carts
        if is_projection_in_triangle(p, vertex[q].carts, u, v):
            face = [q, i, j]
            break
    return face
        
        
        
        
# Is point p in the triangular face?
def is_projection_in_triangle(p, q, u, v):
    # Calculate barycentric coordinates of the projection of p
    n = np.cross(u, v)
    w = p - q
    inv_denom = 1.0 / np.dot(n, n)
    b2 = np.dot( np.cross(u, w), n ) * inv_denom
    b1 = np.dot( np.cross(w, v), n ) * inv_denom
    # Is point in triangular face?
    if b1 >= 0.0 and b2 >= 0.0 and b1+b2 <= 1.0:
        return True
    else:
        return False

    
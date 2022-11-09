# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:38:29 2022

@author: Callum Marples

Precalculate geometric quantities relating to a vertex array, to speed up the
fast marching method.
"""

from leod.fmm_vertex import get_distance
from leod.fmm_vertex import get_angle
from math import sqrt

# Precalculate all neighbour-to-neighbour distances as well as cos(alpha) for
# each face.
def precalculate_grid(vertex):
    n = len(vertex)
    #for i, vi in enumerate(vertex):
    for i in range(n):
        for j in vertex[i].neighbour:
            if i > j:
                # Calculate distance between vertices i and j
                #dist = get_distance(vertex, i, j)

                vertex[i].neighbour[j].distance = sqrt( (vertex[i].carts[0]-vertex[j].carts[0])**2.0 +
                                                        (vertex[i].carts[1]-vertex[j].carts[1])**2.0 +
                                                        (vertex[i].carts[2]-vertex[j].carts[2])**2.0 )
                vertex[j].neighbour[i].distance = vertex[i].neighbour[j].distance
                
    for i in range(n):        
        for j in vertex[i].neighbour:
            for k in vertex[i].neighbour[j].face_angle.keys():
                #cos_alpha = get_angle(vertex, i, j, k)
                
                if vertex[i].neighbour[j].face_angle[k] == 2.0:
                    w1 = vertex[j].carts - vertex[i].carts
                    w2 = vertex[k].carts - vertex[i].carts
                    a = vertex[i].neighbour[j].distance
                    b = vertex[i].neighbour[k].distance
                    vertex[i].neighbour[j].face_angle[k] = (w1[0]*w2[0] + w1[1]*w2[1] + w1[2]*w2[2]) / (a*b)
                    vertex[i].neighbour[k].face_angle[j] = vertex[i].neighbour[j].face_angle[k]
                    
                
                
# Determine all possible faces for each vertex and calculate dot products of
# the edge vectors (i.e. the cosines of the vertex angles)
#def precalcualte_face_angles(vertex):
    #n = len(vertex)
    #for i in range(n):
        
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:38:29 2022

@author: Callum Marples

Precalculate geometric quantities relating to a vertex array, to speed up the
fast marching method.
"""

from leod.fmm_vertex import get_distance
from leod.fmm_vertex import get_angle

# Precalculate all neighbour-to-neighbour distances as well as cos(alpha) for
# each face.
def precalculate_grid(vertex):
    n = len(vertex)
    #for i, vi in enumerate(vertex):
    for i in range(n):
        for j in vertex[i].neighbour:
            if i > j:
                dist = get_distance(vertex, i, j)
            for k in vertex[i].neighbour[j].face_angle.keys():
                cos_alpha = get_angle(vertex, i, j, k)
                
                
                
# Determine all possible faces for each vertex and calculate dot products of
# the edge vectors (i.e. the cosines of the vertex angles)
#def precalcualte_face_angles(vertex):
    #n = len(vertex)
    #for i in range(n):
        
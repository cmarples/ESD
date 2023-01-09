# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:38:29 2022

@author: Callum Marples

Precalculate geometric quantities relating to a vertex array, to speed up the
fast marching method.
"""

import numpy as np
from math import sqrt
from math import acos
from math import pi
from .classes import FmmNeighbour
        
# Precalculate all neighbour-to-neighbour distances as well as cos(alpha) for
# each face.
def grid_precalculate(vertex):
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



def split_obtuse_angles(vertex):
    for i in range(len(vertex)):
        # Find faces to loop over
        face_set = set()
        for j in vertex[i].neighbour.keys():
            for k in vertex[i].neighbour[j].face:
                if k > j:
                    face_set.add((j, k))
        # Loop over each face, splitting if necessary
        for face in face_set:
            j = face[0]
            k = face[1]
            # Check the cosine of the angle to see if it is obtuse.
            if vertex[i].neighbour[j].face_angle[k] < 0.0:
                # Perform split on obtuse angle.
                # Find vertex s to perform split with.
                temp_set = vertex[j].neighbour.keys() & vertex[k].neighbour.keys()
                temp_set.remove(i)
                (s,) = temp_set
                # Add new neighbour
                vertex[i].neighbour[s] = FmmNeighbour()
                # Add outgoing neighbour
                vertex[s].neighbour_set.add(i)
                # Calculate distance
                vertex[i].neighbour[s].distance = sqrt( (vertex[i].carts[0]-vertex[s].carts[0])**2.0 +
                                                        (vertex[i].carts[1]-vertex[s].carts[1])**2.0 +
                                                        (vertex[i].carts[2]-vertex[s].carts[2])**2.0 )
                # Update face completion lists
                for v in range(2):
                    if vertex[i].neighbour[j].face[v] == k:
                        vertex[i].neighbour[j].face[v] = s
                    if vertex[i].neighbour[k].face[v] == j:
                        vertex[i].neighbour[k].face[v] = s
                vertex[i].neighbour[s].face = [j, k]
                # Remove old face angles        
                del vertex[i].neighbour[j].face_angle[k]
                del vertex[i].neighbour[k].face_angle[j]
                # Calculate new face angles
                wj = vertex[j].carts - vertex[i].carts
                wk = vertex[k].carts - vertex[i].carts
                ws = vertex[s].carts - vertex[i].carts
                dj = vertex[i].neighbour[j].distance
                dk = vertex[i].neighbour[k].distance
                ds = vertex[i].neighbour[s].distance
                vertex[i].neighbour[j].face_angle[s] = (wj[0]*ws[0] + wj[1]*ws[1] + wj[2]*ws[2]) / (dj*ds)
                vertex[i].neighbour[s].face_angle[j] = vertex[i].neighbour[j].face_angle[s]
                vertex[i].neighbour[k].face_angle[s] = (wk[0]*ws[0] + wk[1]*ws[1] + wk[2]*ws[2]) / (dk*ds)
                vertex[i].neighbour[s].face_angle[k] = vertex[i].neighbour[k].face_angle[s]
                
                        
# Find number of obtuse angles and the largest angle in the triangulation
def find_obtuse_angles(vertex):
    no_obtuse = 0
    max_angle = 1.0 # Initialise cos(angle) at the largest possible value.
    for i in range(len(vertex)):
        for j in vertex[i].neighbour.keys():
            for k in vertex[i].neighbour[j].face_angle.keys():
                if k > j:
                    if vertex[i].neighbour[j].face_angle[k] < 0.0:
                        no_obtuse += 1
                    if vertex[i].neighbour[j].face_angle[k] < max_angle:
                        max_angle = vertex[i].neighbour[j].face_angle[k]
    max_angle = acos(max_angle) * 180.0 / pi
    return no_obtuse, max_angle
                        
                        
                        
                        
                        
                        
                        
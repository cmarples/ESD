# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 11:03:20 2022

@author: Callum Marples

Generate sphere triangulation using the vertices of an icosahedron.
"""

import math
import numpy as np
import copy
from .fmm_vertex import FmmVertex

# Sphere triangulation using a geodesic polyhedron
def triangulate_sphere(radius=1.0, n=10):

    # Icosahedron coordinates
    t = 0.5 * (1.0 + math.sqrt(5.0)) # Golden ratio
    
    ico_vertex = [ [0, 1, t], [0, -1, t], [0, 1, -t], [0, -1, -t],
                   [1, t, 0], [-1, t, 0], [1, -t, 0], [-1, -t, 0],
                   [t, 0, 1], [-t, 0, 1], [t, 0, -1], [-t, 0, -1] ]
    ico_vertex = np.array(ico_vertex)
    
    # Determine neighbouring vertices (edge length = 2)
    neigh = []
    for i in range(12):
        neigh.append([-1] * 5)
        neighbour_no = 0
        for j in range(12):
            if i != j:
                d2 = (ico_vertex[i][0]-ico_vertex[j][0])**2 + (ico_vertex[i][1]-ico_vertex[j][1])**2 + (ico_vertex[i][2]-ico_vertex[j][2])**2
                if math.fabs(d2 - 4.0) < 1e-9:
                    neigh[i][neighbour_no] = j
                    neighbour_no += 1
    
    # Determine faces by considering each trio in turn
    ico_face = []
    for i in range(10):
        for j in range(i, 11):
            for k in range(j, 12):
                if j in neigh[i] and k in neigh[i] and i in neigh[j] and k in neigh[j] and i in neigh[k] and j in neigh[k]:
                    ico_face.append([i, j, k])
    
    vertex_carts = []
    face = []
    
    # Determine triangle vertices and faces
    for i in range(20):
        # Get icosahedron vertices
        v0 = ico_vertex[ico_face[i][0]]
        v1 = ico_vertex[ico_face[i][1]]
        v2 = ico_vertex[ico_face[i][2]]
        # Triangle vertices
        count = 0
        v012 = []
        for j in range(n+1):
            v01 = (v0*(n-j) + v1*(j)) / n
            v02 = (v0*(n-j) + v2*(j)) / n
            for k in range(j+1):
                count += 1
                if j == 0:
                    v012.append(v01)
                else:
                    v012.append((v01*(j-k) + v02*k) / j)
                
        # Triangle faces
        count = 0
        face_tri = []
        for j in range(1, n+1):
            for k in range(1, j+1):
                count += 1
                no = int(j*(j-1)/2 + k) - 1
                face_tri.append([no, no+j, no+j+1])
        for j in range(2, n+1):
            for k in range(1, j):
                count += 1
                no = int(j*(j-1)/2 + k) - 1
                face_tri.append([no, no+j+1, no+1])
                
        # Merge faces and vertices
        m = len(vertex_carts)
        for j in range(len(v012)):
            vertex_carts.append(v012[j])
        for j in range(len(face_tri)):
            face_temp = [face_tri[j][0]+m, face_tri[j][1]+m, face_tri[j][2]+m]
            face.append(face_temp)
    
    # Remove duplicates
    q = len(vertex_carts)
    q2 = copy.copy(q)
    for i in range(q-1):
        jj = i
        for j in range(i+1, q+1):
            jj += 1
            #print([i, j])
            #if np.linalg.norm(vertex_carts[i]-vertex_carts[jj]) < 1e-5:
            if ( (math.fabs(vertex_carts[i][0]-vertex_carts[jj][0]) < 1e-5) and
                 (math.fabs(vertex_carts[i][1]-vertex_carts[jj][1]) < 1e-5) and
                 (math.fabs(vertex_carts[i][2]-vertex_carts[jj][2]) < 1e-5) ):
                del vertex_carts[jj]

                q2 -= 1
                for k in range(len(face)):
                    for l in range(3):
                        if face[k][l] == jj:
                            face[k][l] = i
                        elif face[k][l] > jj:
                            face[k][l] = face[k][l] - 1
                jj -= 1
            if jj >= q2-1:
                break
        if i >= q2-2:
            break
    
    # Project onto sphere surface
    for i in range(q2):
        vertex_carts[i] *= radius / np.linalg.norm(vertex_carts[i])
                
    # Create list of FmmVertex.
    vertex = []
    for i in range(q2):
        vertex.append(FmmVertex(i, vertex_carts[i]))
        # Find neighbours using the faces.
        edge_list = []
        for j in range(len(face)):
            if i in face[j]:
                for k in range(3):
                    edge_list.append(face[j][k])
        for j in edge_list:
            if (j not in vertex[i].neighbour) and j != i:
                vertex[i].neighbour.append(j)
                
    # Get faces belonging to each vertex
    for i in range(len(face)):
        vertex[face[i][0]].face.append([face[i][1], face[i][2]])
        vertex[face[i][1]].face.append([face[i][0], face[i][2]])
        vertex[face[i][2]].face.append([face[i][0], face[i][1]])
        
    return vertex
        
        
   
        

    
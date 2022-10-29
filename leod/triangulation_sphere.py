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
    
    # Determine faces by considering each trio in turn.
    ico_face = []
    for i in range(10):
        for j in range(i, 11):
            for k in range(j, 12):
                if j in neigh[i] and k in neigh[i] and i in neigh[j] and k in neigh[j] and i in neigh[k] and j in neigh[k]:
                    ico_face.append([i, j, k])
    
    # Determine the five faces surrounding each vertex.
    ico_vertex_face = []
    for j in range(12):
        ico_vertex_face.append([])
        for i in range(20):
            if j in ico_face[i]:
                ico_vertex_face[j].append(i)
                
    # Determine the two faces either side of each edge.
    ico_edge_face = {}
    for i in range(12):
        for j in neigh[i]:
            if j > i:
                ico_edge_face[(i, j)] = []
                for k in range(20):
                    if i in ico_face[k]:
                        if j in ico_face[k]:
                            ico_edge_face[(i, j)].append(k)
    
    
    
    def is_corner_vertex_defined(face_no, vertex_no):
        is_defined = False
        for f in ico_vertex_face[vertex_no]:
            if f < face_no:
                is_defined = True
                break
        return is_defined
    
    def is_edge_vertex_defined(face_no, vi, vj):
        is_defined = False
        if vi < vj:
            x = (vi, vj)
        else:
            x = (vj, vi)
        f1, f2 = ico_edge_face[x]
        if f1 < face_no or f2 < face_no:
            is_defined = True
        return is_defined
    
    def process_corner_vertex(vertex_no, index):
        if ico_index[vertex_no] == -1:
            vertex.append(FmmVertex(index, ico_vertex[vertex_no]))
            vertex_map[(j, k)] = index
            ico_index[vertex_no] = index
            index += 1
        return index
    
    def process_edge_vertex(vx, vy, jk, index):
        if vx < vy:
            x_lt_y = True
            val = ico_edge_index[(vx, vy)][jk-1]
        else:
            x_lt_y = False
            val = ico_edge_index[(vx, vy)][n-1-jk]                   
        if val == -1:
            vertex.append(FmmVertex(index, (v01*(j-k) + v02*k) / j))
            vertex_map[(j, k)] = index
            if x_lt_y == True:
                ico_edge_index[(vx, vy)][jk-1] = index
            else:
                ico_edge_index[(vx, vy)][n-1-jk] = index
            index += 1 
        return index
                                
                                
    vertex = []
    ico_index = [-1] * 12 # Map from icosahedral index to triangulation index
    ico_edge_index = {}   # Map from an edge of a icosahedral face to an index
    for i in range(12):
        for j in neigh[i]:
            if j > i:
                ico_edge_index[(i, j)] = [-1] * (n-1)
    index = 0
    
    
            
    # Determine vertices, faces and neighbours of the triangulation
    for i in range(20):
        
        # Get icosahedron vertices
        v0 = ico_vertex[ico_face[i][0]]
        v1 = ico_vertex[ico_face[i][1]]
        v2 = ico_vertex[ico_face[i][2]]

        # Triangle vertices
        vertex_map = {}
        
        for j in range(n+1):
            
            v01 = (v0*(n-j) + v1*(j)) / n
            v02 = (v0*(n-j) + v2*(j)) / n
            
            if j == 0:
                # Vertex v0 may have already been defined.
                index = process_corner_vertex(ico_face[i][0], index)
            else:
                for k in range(j+1):
                    if j == n:
                        if k == 0:
                            # Vertex v1 may have already been defined.
                            index = process_corner_vertex(ico_face[i][1], index) 
                        elif k == j:
                            # Vertex v2 may have already been defined.
                            index = process_corner_vertex(ico_face[i][2], index)
                        else:
                            # Vertex along v1-v2 line may have already been defined
                            index = process_edge_vertex(ico_face[i][1], ico_face[i][2], k, index)
                    else:
                        if k == 0:
                            # Vertex along v0-v1 line may have already been defined
                            index = process_edge_vertex(ico_face[i][0], ico_face[i][1], j, index)    
                        elif k == j:
                            # Vertex along v0-v2 line may have already been defined
                            index = process_edge_vertex(ico_face[i][0], ico_face[i][2], j, index)
                        else:
                            # Vertex is within the icosahedral face
                            vertex.append(FmmVertex(index, (v01*(j-k) + v02*k) / j))
                            vertex_map[(j, k)] = index
                            index += 1
                
                       
        # Find neighbours and faces for each vertex
        x=0
        
        
        
        
        '''        
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
    '''
    # Project onto sphere surface
    for i in range(len(vertex)):
        vertex[i].carts *= radius / np.linalg.norm(vertex[i].carts)
        
    return vertex
'''            
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
            if (j not in vertex[i].distance_to_neighbour.keys()) and j != i:
                vertex[i].distance_to_neighbour[j] = -1.0
                
    # Get faces belonging to each vertex
    for i in range(len(face)):
        vertex[face[i][0]].face.append([face[i][1], face[i][2]])
        vertex[face[i][1]].face.append([face[i][0], face[i][2]])
        vertex[face[i][2]].face.append([face[i][0], face[i][1]])
'''        
    
        
        
   
        

    
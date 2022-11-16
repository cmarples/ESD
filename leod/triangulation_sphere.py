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
from .fmm_vertex import FmmNeighbour

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
    
    def process_corner_vertex(vertex_no, index):
        if ico_index[vertex_no] == -1:
            vertex.append(FmmVertex(index, ico_vertex[vertex_no]))
            vertex[-1].face = []
            if j == 0:
                vertex_map[(j, 0)] = index
            else:
                vertex_map[(j, k)] = index
            ico_index[vertex_no] = index
            index += 1
        else:
            if j == 0:
                vertex_map[(j, 0)] = ico_index[vertex_no]
            else:
                vertex_map[(j, k)] = ico_index[vertex_no]
        return index
    
    def process_edge_vertex(vx, vy, jk, index):
        if vx < vy:
            x_lt_y = True
            vertex_index = ico_edge_index[(vx, vy)][jk-1]
        else:
            x_lt_y = False
            vertex_index = ico_edge_index[(vx, vy)][n-1-jk]                   
        if vertex_index == -1:
            vertex.append(FmmVertex(index, (v01*(j-k) + v02*k) / j))
            vertex[-1].face = []
            vertex_map[(j, k)] = index
            if x_lt_y == True:
                ico_edge_index[(vx, vy)][jk-1] = index
            else:
                ico_edge_index[(vx, vy)][n-1-jk] = index
            index += 1 
        else:
            vertex_map[(j, k)] = vertex_index
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
            
            if j == 0: # Vertex v0 may have already been defined.
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
                            vertex[-1].face = []
                            vertex_map[(j, k)] = index
                            index += 1
                
                       
        # Find neighbours and faces for each vertex
        for j in range(n+1):
            if j == 0: # Vertex v0
                v = vertex_map[(j, 0)]
                vertex[v].neighbour[vertex_map[(1, 0)]] = FmmNeighbour()
                vertex[v].neighbour[vertex_map[(1, 1)]] = FmmNeighbour()
                vertex[v].face.append([vertex_map[(1, 1)], vertex_map[(1, 0)]])
            else:
                for k in range(j+1):
                    v = vertex_map[(j, k)]
                    if j == n:
                        if k == 0: # Vertex v1 
                            vertex[v].neighbour[vertex_map[(n-1, 0)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(n, 1)]] = FmmNeighbour()
                            vertex[v].face.append([vertex_map[(n-1, 0)], vertex_map[(n, 1)]])
                        elif k == j: # Vertex v2
                            vertex[v].neighbour[vertex_map[(n-1, n-1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(n, n-1)]] = FmmNeighbour()
                            vertex[v].face.append([vertex_map[(n-1, n-1)], vertex_map[(n, n-1)]])
                        else:
                            # Vertex along v1-v2 line
                            vertex[v].neighbour[vertex_map[(n-1, k)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(n, k+1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(n, k-1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(n-1, k-1)]] = FmmNeighbour()
                            vertex[v].face.append([vertex_map[(n-1, k)], vertex_map[(n, k+1)]])
                            vertex[v].face.append([vertex_map[(n, k-1)], vertex_map[(n-1, k-1)]])
                            vertex[v].face.append([vertex_map[(n-1, k-1)], vertex_map[(n-1, k)]])
                    else:
                        if k == 0:
                            # Vertex along v0-v1 line
                            vertex[v].neighbour[vertex_map[(j-1, 0)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j, 1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j+1, 1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j+1, 0)]] = FmmNeighbour()
                            vertex[v].face.append([vertex_map[(j-1, 0)], vertex_map[(j, 1)]])
                            vertex[v].face.append([vertex_map[(j, 1)], vertex_map[(j+1, 1)]])
                            vertex[v].face.append([vertex_map[(j+1, 1)], vertex_map[(j+1, 0)]])
                        elif k == j:
                            # Vertex along v0-v2 line
                            vertex[v].neighbour[vertex_map[(j+1, j+1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j+1, j)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j, j-1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j-1, j-1)]] = FmmNeighbour()
                            vertex[v].face.append([vertex_map[(j+1, j+1)], vertex_map[(j+1, j)]])
                            vertex[v].face.append([vertex_map[(j+1, j)], vertex_map[(j, j-1)]])
                            vertex[v].face.append([vertex_map[(j, j-1)], vertex_map[(j-1, j-1)]])
                        else:
                            # Vertex is within the icosahedral face
                            vertex[v].neighbour[vertex_map[(j-1, k)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j, k+1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j+1, k+1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j+1, k)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j, k-1)]] = FmmNeighbour()
                            vertex[v].neighbour[vertex_map[(j-1, k-1)]] = FmmNeighbour()
                            vertex[v].face.append([vertex_map[(j-1, k)], vertex_map[(j, k+1)]])
                            vertex[v].face.append([vertex_map[(j, k+1)], vertex_map[(j+1, k+1)]])
                            vertex[v].face.append([vertex_map[(j+1, k+1)], vertex_map[(j+1, k)]])
                            vertex[v].face.append([vertex_map[(j+1, k)], vertex_map[(j, k-1)]])
                            vertex[v].face.append([vertex_map[(j, k-1)], vertex_map[(j-1, k-1)]])
                            vertex[v].face.append([vertex_map[(j-1, k-1)], vertex_map[(j-1, k)]])
    
    # Find possible faces for each neighbouring pair
    cos_alpha = 2.0 # Initialisation for face angles
    for i in range(len(vertex)):
        for j in vertex[i].neighbour.keys():
            k_list = []
            for face_i in vertex[i].face:
                if j in face_i:
                    if face_i[0] == j:
                        k_list.append(face_i[1])
                    else:
                        k_list.append(face_i[0])
            vertex[i].neighbour[j].face = k_list
            vertex[i].neighbour[j].face_angle[k_list[0]] = cos_alpha
            vertex[i].neighbour[j].face_angle[k_list[1]] = cos_alpha
                    
                            
                            
                            
    # Project onto sphere surface
    for i in range(len(vertex)):
        vertex[i].carts *= radius / np.linalg.norm(vertex[i].carts)
        
    return vertex

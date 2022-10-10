# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:49:51 2022

@author: Callum Marples

Generate FmmGrid using an ellipsoid with scaled spherical polar coordinates
"""

import math
import numpy as np

from .ellipsoid_shape import EllipsoidShape
from .fmm_vertex import FmmVertex

def generate_polar_graph(shape, no_theta, no_phi, is_connect_8=False, is_Dijkstra=False):
    
    # Create list of FmmVertex
    vertex = []
    
    # Structured grid information
    no_vertices = (no_theta - 2)*no_phi + 2
    delta_theta = math.pi / (no_theta - 1)
    delta_phi = 2.0*math.pi / no_phi
    
    # Compute lists of theta and phi values
    theta_list = [0.0] * no_theta
    phi_list = [0.0] * (no_phi + 1)
    for i in range(no_theta):
        theta_list[i] = i *delta_theta
    for i in range(no_phi+1):
        phi_list[i] = i * delta_phi
    
    # Construct the vertex array in the FmmGrid
    polar_index = []
    for i in range(no_vertices):
        
        # theta and phi indices
        th_index = get_theta_index(i, no_vertices, no_theta, no_phi)
        ph_index = get_phi_index(i, th_index, no_phi)
        polar_index.append([th_index, ph_index])
        
        # Cartesian coordinates
        carts = np.array( shape.polar2cart(theta_list[th_index], phi_list[ph_index]) )
        
        # Create new vertex
        vertex.append(FmmVertex(i, carts))
        
    # Find neighbouring vertices and update the graph
    for i in range(no_vertices):    
        find_neighbour_indices(vertex, i, polar_index[i][0], polar_index[i][1], 
                               no_theta, no_phi, no_vertices, is_connect_8)
    
    # Find Fast Marching update triangles
    if is_Dijkstra == False:
        for i in range(no_vertices):
            find_update_triangles(vertex, i, no_theta, no_phi, no_vertices, is_connect_8)
    
    
    
    
    
    
    return vertex

############################# Define subroutines #############################

# Get vertex index from theta and phi indices
def get_vertex_index(theta_index, phi_index, no_vertices, no_theta, no_phi):
    if theta_index > 0 and theta_index < no_theta-1:
        return 1 + phi_index + no_phi*(theta_index-1)
    elif theta_index == 0:
        return 0
    else: # theta_index = no_theta - 1
        return no_vertices - 1
        
# Get theta index from vertex index
def get_theta_index(vertex_index, no_vertices, no_theta, no_phi):
    if vertex_index > 0 and vertex_index < no_vertices-1:
        return math.ceil(float(vertex_index)/float(no_phi))
    elif vertex_index == 0:
        return 0
    elif vertex_index == no_vertices-1:
        return no_theta-1
    else:
        return no_theta
        
# Get phi index from vertex index and theta index
def get_phi_index(vertex_index, theta_index, no_phi):
    return ( vertex_index - 1 - no_phi*(theta_index-1) )

# Find indices of neighbouring vertices
# Each pixel has four neighbours. For non-polar vertices, neighbours are
# always ordered as "Up, Down, Left, Right".
def find_neighbour_indices(vertex, i, th, ph, no_theta, no_phi, no_vertices, is_connect_8=False):
    
    # Pole neighbours
    if i == 0:
        
        # North pole                  
        for k in range(no_phi):
            vertex[0].neighbour.append(1+k)
            
    elif i == no_vertices-1: 
        
        # South pole
        for k in range(no_phi):
            #vertex[no_vertices-1].neighbour.append(no_vertices-2-k)
            vertex[no_vertices-1].neighbour.append(no_vertices-1-no_phi+k)
            
    else:     
        # Not a pole

        # Up neighbour (minus theta)
        if th == 1:
            vertex_neighbour = 0
        else:
            vertex_neighbour = get_vertex_index(th-1, ph, no_vertices, no_theta, no_phi)
        vertex[i].neighbour.append(vertex_neighbour)
        
        # Down neighbour (plus theta)
        if th == no_theta-2:
            vertex_neighbour = no_vertices - 1
        else:
            vertex_neighbour = get_vertex_index(th+1, ph, no_vertices, no_theta, no_phi)
        vertex[i].neighbour.append(vertex_neighbour) 
        
        # Left neighbour (minus phi)
        if ph == 0:
            vertex_neighbour = get_vertex_index(th, no_phi-1, no_vertices, no_theta, no_phi)
        else:
            vertex_neighbour = get_vertex_index(th, ph-1, no_vertices, no_theta, no_phi)
        vertex[i].neighbour.append(vertex_neighbour)
        
        # Right neighbour (plus phi)
        if ph == no_phi-1:
            vertex_neighbour = get_vertex_index(th, 0, no_vertices, no_theta, no_phi)
        else:
            vertex_neighbour = get_vertex_index(th, ph+1, no_vertices, no_theta, no_phi)
        vertex[i].neighbour.append(vertex_neighbour) 
        
        # Diagonal neighbours
        if is_connect_8 == True:
            # Up-left neighbour (minus theta, minus phi)
            if th == 1:
                vertex_neighbour = 0
            else:
                if ph == 0:
                    vertex_neighbour = get_vertex_index(th-1, no_phi-1, no_vertices, no_theta, no_phi)
                else:
                    vertex_neighbour = get_vertex_index(th-1, ph-1, no_vertices, no_theta, no_phi)
            vertex[i].neighbour.append(vertex_neighbour)
            # Up-right neighbour (minus theta, plus phi)
            if th == 1:
                vertex_neighbour = 0
            else:
                if ph == no_phi-1:
                    vertex_neighbour = get_vertex_index(th-1, 0, no_vertices, no_theta, no_phi)
                else:
                    vertex_neighbour = get_vertex_index(th-1, ph+1, no_vertices, no_theta, no_phi)
            vertex[i].neighbour.append(vertex_neighbour)
            # Down-left neighbour (plus theta, minus phi)
            if th == no_theta-2:
                vertex_neighbour = no_vertices - 1
            else:
                if ph == 0:
                    vertex_neighbour = get_vertex_index(th+1, no_phi-1, no_vertices, no_theta, no_phi)
                else:
                    vertex_neighbour = vertex_neighbour = get_vertex_index(th+1, ph-1, no_vertices, no_theta, no_phi)
            vertex[i].neighbour.append(vertex_neighbour)
            # Down-right neighbour (plus theta, plus phi)
            if th == no_theta-2:
                vertex_neighbour = no_vertices - 1
            else:
                if ph == no_phi-1:
                    vertex_neighbour = get_vertex_index(th+1, 0, no_vertices, no_theta, no_phi)
                else:
                    vertex_neighbour = get_vertex_index(th+1, ph+1, no_vertices, no_theta, no_phi)
            vertex[i].neighbour.append(vertex_neighbour)

# Get the two adjacent phi indices (accounting for the phi boundary)
def get_adjacent_phi(j, no_phi):
    j_plus = j + 1
    j_minus = j - 1
    if j_plus == no_phi:
        j_plus -= no_phi
    if j_minus == -1:
        j_minus += no_phi
    return [j_plus, j_minus]

# Find the possible triangles that can be used to update vertex i
def find_update_triangles(vertex, i, no_theta, no_phi, no_vertices, is_connect_8):
    
    # Handle poles seperately
    if i == 0:
        # North pole
        for j in range(no_phi-1):
            vertex[i].face.append([j+1, j+2])
        vertex[i].face.append([no_phi, 1])
    elif i == no_vertices-1:
        # South pole
        k = no_vertices-1-no_phi # Starting vertex
        for j in range(no_phi-1):
            
            #[j_plus, j_minus] = get_adjacent_phi(j, no_phi)
            #vertex[i].neighbour_faces.append([no_vertices-1-no_phi+j_plus, no_vertices-1-no_phi+j_minus])
            vertex[i].face.append([k+j, k+j+1])
        vertex[i].face.append([no_vertices-2, k])
    else:    
        if is_connect_8 == False:
            # Non-polar vertex with 4 neighbours
            vertex[i].face.append([vertex[i].neighbour[0], vertex[i].neighbour[3]]) # Up-right
            vertex[i].face.append([vertex[i].neighbour[3], vertex[i].neighbour[1]]) # Down-right
            vertex[i].face.append([vertex[i].neighbour[1], vertex[i].neighbour[2]]) # Down-left
            vertex[i].face.append([vertex[i].neighbour[2], vertex[i].neighbour[0]]) # Up-left
            
        else:
            th_index = get_theta_index(i, no_vertices, no_theta, no_phi)
            
            if th_index == 1:
                # North pole adjacent
                del vertex[i].neighbour[4]
                del vertex[i].neighbour[4] # Remove duplicated north pole neighbours
                vertex[i].face.append([vertex[i].neighbour[0], vertex[i].neighbour[3]]) # Up-right
                vertex[i].face.append([vertex[i].neighbour[3], vertex[i].neighbour[5]]) # Down-right-right
                vertex[i].face.append([vertex[i].neighbour[5], vertex[i].neighbour[1]]) # Down-right-left
                vertex[i].face.append([vertex[i].neighbour[1], vertex[i].neighbour[4]]) # Down-left-right
                vertex[i].face.append([vertex[i].neighbour[4], vertex[i].neighbour[2]]) # Down-left-left
                vertex[i].face.append([vertex[i].neighbour[2], vertex[i].neighbour[0]]) # Up-left
                
            elif th_index == no_theta - 2:
                # South pole adjacent
                del vertex[i].neighbour[6]
                del vertex[i].neighbour[6] # Remove duplicated south pole neighbours
                vertex[i].face.append([vertex[i].neighbour[0], vertex[i].neighbour[5]]) # Up-right-left
                vertex[i].face.append([vertex[i].neighbour[5], vertex[i].neighbour[3]]) # Up-right-right
                vertex[i].face.append([vertex[i].neighbour[3], vertex[i].neighbour[1]]) # Down-right
                vertex[i].face.append([vertex[i].neighbour[1], vertex[i].neighbour[2]]) # Down-left
                vertex[i].face.append([vertex[i].neighbour[2], vertex[i].neighbour[4]]) # Up-left-left
                vertex[i].face.append([vertex[i].neighbour[4], vertex[i].neighbour[0]]) # Up-left-right
                
            else:
                
                # Non-polar vertex with 8 neighbours
                vertex[i].face.append([vertex[i].neighbour[0], vertex[i].neighbour[5]]) # Up-right-left
                vertex[i].face.append([vertex[i].neighbour[5], vertex[i].neighbour[3]]) # Up-right-right
                vertex[i].face.append([vertex[i].neighbour[3], vertex[i].neighbour[7]]) # Down-right-right
                vertex[i].face.append([vertex[i].neighbour[7], vertex[i].neighbour[1]]) # Down-right-left
                vertex[i].face.append([vertex[i].neighbour[1], vertex[i].neighbour[6]]) # Down-left-right
                vertex[i].face.append([vertex[i].neighbour[6], vertex[i].neighbour[2]]) # Down-left-left
                vertex[i].face.append([vertex[i].neighbour[2], vertex[i].neighbour[4]]) # Up-left-left
                vertex[i].face.append([vertex[i].neighbour[4], vertex[i].neighbour[0]]) # Up-left-right

# Check that each triangle is valid by calculating each vertex angle
def check_triangles(vertex):
    no_obtuse = 0
    # Check each vertex in turn
    for i in range(len(vertex)):
        vertex[i].face_valid = []
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
            if cos_alpha > 0 or math.fabs(cos_alpha) < 1.0e-15:
                vertex[i].face_valid.append(True)
            else:
                vertex[i].face_valid.append(False) 
                no_obtuse += 1
    return no_obtuse
        

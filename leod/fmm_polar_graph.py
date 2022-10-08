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
    
    # Check that these triangles are acute (and thus valid for the FMM)
    #check_triangles(vertex)
    
    
    
    
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
            vertex[0].neighbour_distance.append(-1.0)
            
    elif i == no_vertices-1: 
        
        # South pole
        for k in range(no_phi):
            vertex[no_vertices-1].neighbour.append(no_vertices-2-k)
            vertex[no_vertices-1].neighbour_distance.append(-1.0)
            
    else:     
        # Not a pole
                         
        # Initialise distances
        if is_connect_8 == True:
            vertex[i].neighbour_distance = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
        else:
            vertex[i].neighbour_distance = [-1.0, -1.0, -1.0, -1.0]
            
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
        for j in range(no_phi):
            [j_plus, j_minus] = get_adjacent_phi(j, no_phi)
            vertex[i].neighbour_faces.append([j_plus+1, j_minus+1])
    elif i == no_vertices-1:
        # South pole
        for j in range(no_phi):
            [j_plus, j_minus] = get_adjacent_phi(j, no_phi)
            vertex[i].neighbour_faces.append([no_vertices-1-no_phi+j_plus, no_vertices-1-no_phi+j_minus])
    else:    
        if is_connect_8 == False:
            # Non-polar vertex with 4 neighbours
            
            # Up     (need left and right)
            vertex[i].neighbour_faces.append([vertex[i].neighbour[2], vertex[i].neighbour[3]])
            # Down   (need left and right)
            vertex[i].neighbour_faces.append([vertex[i].neighbour[2], vertex[i].neighbour[3]])
            # Left   (need up and down)
            vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[1]])
            # Right  (need up and down)
            vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[1]])
                    
        else:
            th_index = get_theta_index(i, no_vertices, no_theta, no_phi)
            # Non-polar vertex with 8 neighbours
            if th_index == 1:
                # North pole adjacent
                del vertex[i].neighbour[4]
                del vertex[i].neighbour[4] # Remove duplicated north pole neighbours
                # Up         (need left and right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[2], vertex[i].neighbour[3]])
                # Down       (need down-left and down-right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[4], vertex[i].neighbour[5]])
                # Left       (need up and down-left)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[4]])
                # Right      (need up and down-right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[5]])
                # Down-left  (need down and left)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[1], vertex[i].neighbour[2]])
                # Down-right (need down and right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[1], vertex[i].neighbour[3]])
                
            elif th_index == no_theta - 2:
                # South pole adjacent
                del vertex[i].neighbour[6]
                del vertex[i].neighbour[6] # Remove duplicated south pole neighbours
                # Up         (need up-left and up-right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[4], vertex[i].neighbour[5]])
                # Down       (need left and right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[2], vertex[i].neighbour[3]])
                # Left       (need down and up-left)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[1], vertex[i].neighbour[4]])
                # Right      (need down and up-right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[1], vertex[i].neighbour[5]])
                # Up-left    (need up and left)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[2]])
                # Up-right   (need up and right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[3]])
                
            else:
                # Up         (need up-left and up-right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[4], vertex[i].neighbour[5]])
                # Down       (need down-left and down-right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[6], vertex[i].neighbour[7]])
                # Left       (need up-left and down-left)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[4], vertex[i].neighbour[6]])
                # Right      (need up-right and down-right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[5], vertex[i].neighbour[7]])
                # Up-left    (need up and left)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[2]])
                # Up-right   (need up and right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[0], vertex[i].neighbour[3]])
                # Down-left  (need down and left)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[1], vertex[i].neighbour[2]])
                # Down-right (need down and right)
                vertex[i].neighbour_faces.append([vertex[i].neighbour[1], vertex[i].neighbour[3]])
            
# Check that each triangle is valid by calculating each vertex angle
def check_triangles(vertex):
    
    # Check each vertex in turn
    for i in range(len(vertex)):
        vertex[i].face_valid = []
        vertex[i].face_angle = []
        #if i == 0:
        #    x = 1
        #elif i == len(vertex)-1:
        #    x=1
        #else:
        for j_count in range(len(vertex[i].neighbour)):
            j = vertex[i].neighbour[j_count]
            w1 = vertex[i].carts - vertex[j].carts
            vertex[i].face_valid.append([False, False])
            vertex[i].face_angle.append([-2.0, -2.0])
            for k_count in range(2):
                k = vertex[i].neighbour_faces[j_count][k_count]
                w2 = vertex[i].carts - vertex[k].carts
                # Need angle between vectors w1 and w2 to be <= 90 degrees
                cos_alpha = np.dot(w1, w2)
                vertex[i].face_angle[j_count][k_count] = cos_alpha
                if cos_alpha > 0 or math.fabs(cos_alpha) < 1.0e-15:
                    vertex[i].face_valid[j_count][k_count] = True
    
    
    
    
    
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:49:51 2022

@author: Callum Marples

Generate FmmGrid using an ellipsoid with scaled spherical polar coordinates
"""

import math
import numpy as np

from .fmm_vertex import FmmVertex
from .fmm_vertex import FmmNeighbour
from .binary_search import binary_search

class PolarGrid:
    def __init__(self, no_theta, no_phi):
        self.no_theta = no_theta
        self.no_phi = no_phi
        self.no_vertices = (no_theta - 2)*no_phi + 2
        self.delta_theta = math.pi / (no_theta - 1)
        self.delta_phi = 2.0*math.pi / no_phi
        
        # Compute lists of theta and phi values
        self.theta_list = [0.0] * no_theta
        self.phi_list = [0.0] * (no_phi + 1)
        for i in range(no_theta):
            self.theta_list[i] = i * self.delta_theta
        for i in range(no_phi+1):
            self.phi_list[i] = i * self.delta_phi
    

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
    #if is_Dijkstra == False:
    #    for i in range(no_vertices):
    #        find_update_triangles(vertex, i, no_theta, no_phi, no_vertices, is_connect_8)
    
    
    
    
    
    
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

# Find theta index closest to given theta value, th
def find_theta_index(theta_list, th):
    if th < math.pi:
        return binary_search(theta_list, th)
    else:
        return len(theta_list)-1

# Find phi index closest to given phi value, ph
def find_phi_index(phi_list, ph):
    index = binary_search(phi_list, ph)
    if index == len(phi_list)-1:
        return 0
    else:
        return index

# Find index of vertex closest to th and ph
def find_vertex_index(theta_list, phi_list, th, ph):
    th_index = find_theta_index(theta_list, th)
    ph_index = find_phi_index(phi_list, ph)
    no_vertices = (len(theta_list) - 2)*(len(phi_list)-1) + 2
    return [get_vertex_index(th_index, ph_index, no_vertices, len(theta_list), len(phi_list)-1), th_index,ph_index]
    
# Find indices of neighbouring vertices
# Each pixel has four neighbours. For non-polar vertices, neighbours are
# always ordered as "Up, Down, Left, Right".
def find_neighbour_indices(vertex, i, th, ph, no_theta, no_phi, no_vertices, is_connect_8=False):
    
    k0 = no_vertices-1-no_phi # First vertex at south pole adjacent band  
    cos_alpha = 2.0 # Initialisation for face angles
    
    # Pole neighbours
    if i == 0:
        
        # North pole
        vertex[0].neighbour[1] = FmmNeighbour()
        vertex[0].neighbour[1].face = [no_phi, 2]
        vertex[0].neighbour[1].face_angle[no_phi] = cos_alpha
        vertex[0].neighbour[1].face_angle[2] = cos_alpha               
        for k in range(1, no_phi-1):
            #vertex[0].neighbour.append(1+k)
            temp = 1+k
            vertex[0].neighbour[temp] = FmmNeighbour()
            vertex[0].neighbour[temp].face = [k, 2+k]
            vertex[0].neighbour[temp].face_angle[k] = cos_alpha
            vertex[0].neighbour[temp].face_angle[2+k] = cos_alpha 
        vertex[0].neighbour[no_phi] = FmmNeighbour()
        vertex[0].neighbour[no_phi].face = [no_phi-1, 1]
        vertex[0].neighbour[no_phi].face_angle[no_phi-1] = cos_alpha
        vertex[0].neighbour[no_phi].face_angle[1] = cos_alpha
            
    elif i == no_vertices-1: 
        
        # South pole
        vertex[no_vertices-1].neighbour[k0] = FmmNeighbour()
        vertex[no_vertices-1].neighbour[k0].face = [no_vertices-2, k0+1]
        vertex[no_vertices-1].neighbour[k0].face_angle[no_vertices-2] = cos_alpha
        vertex[no_vertices-1].neighbour[k0].face_angle[k0+1] = cos_alpha
        for k in range(1, no_phi-1):
            #vertex[no_vertices-1].neighbour.append(no_vertices-2-k)
            #vertex[no_vertices-1].neighbour.append(no_vertices-1-no_phi+k)
            temp = k0+k
            vertex[no_vertices-1].neighbour[temp] = FmmNeighbour()
            vertex[no_vertices-1].neighbour[temp].face = [temp-1, temp+1]
            vertex[no_vertices-1].neighbour[temp].face_angle[temp-1] = cos_alpha
            vertex[no_vertices-1].neighbour[temp].face_angle[temp+1] = cos_alpha  
        vertex[no_vertices-1].neighbour[no_vertices-2] = FmmNeighbour()
        vertex[no_vertices-1].neighbour[no_vertices-2].face = [no_vertices-3, k0]
        vertex[no_vertices-1].neighbour[no_vertices-2].face_angle[no_vertices-3] = cos_alpha
        vertex[no_vertices-1].neighbour[no_vertices-2].face_angle[k0] = cos_alpha
        
    else:     
        # Not a pole

        # Up neighbour (minus theta)
        if th == 1:
            j_up = 0
        else:
            j_up = get_vertex_index(th-1, ph, no_vertices, no_theta, no_phi)
        #vertex[i].neighbour.append(vertex_neighbour)
        vertex[i].neighbour[j_up] = FmmNeighbour()    
        
        # Down neighbour (plus theta)
        if th == no_theta-2:
            j_dn = no_vertices - 1
        else:
            j_dn = get_vertex_index(th+1, ph, no_vertices, no_theta, no_phi)
        #vertex[i].neighbour.append(vertex_neighbour)
        vertex[i].neighbour[j_dn] = FmmNeighbour() 
        
        # Left neighbour (minus phi)
        if ph == 0:
            j_lt = get_vertex_index(th, no_phi-1, no_vertices, no_theta, no_phi)
        else:
            j_lt = get_vertex_index(th, ph-1, no_vertices, no_theta, no_phi)
        #vertex[i].neighbour.append(vertex_neighbour)
        vertex[i].neighbour[j_lt] = FmmNeighbour()
        
        # Right neighbour (plus phi)
        if ph == no_phi-1:
            j_rt = get_vertex_index(th, 0, no_vertices, no_theta, no_phi)
        else:
            j_rt = get_vertex_index(th, ph+1, no_vertices, no_theta, no_phi)
        #vertex[i].neighbour.append(vertex_neighbour)
        vertex[i].neighbour[j_rt] = FmmNeighbour()
        
        # Diagonal neighbours
        if is_connect_8 == True:
            if i > no_phi:
                # Up-left neighbour (minus theta, minus phi)
                if th == 1:
                    j_uplt = 0
                else:
                    if ph == 0:
                        j_uplt = get_vertex_index(th-1, no_phi-1, no_vertices, no_theta, no_phi)
                    else:
                        j_uplt = get_vertex_index(th-1, ph-1, no_vertices, no_theta, no_phi)
                #vertex[i].neighbour.append(vertex_neighbour)
                vertex[i].neighbour[j_uplt] = FmmNeighbour()
                
                # Up-right neighbour (minus theta, plus phi)
                if th == 1:
                    j_uprt = 0
                else:
                    if ph == no_phi-1:
                        j_uprt = get_vertex_index(th-1, 0, no_vertices, no_theta, no_phi)
                    else:
                        j_uprt = get_vertex_index(th-1, ph+1, no_vertices, no_theta, no_phi)
                #vertex[i].neighbour.append(vertex_neighbour)
                vertex[i].neighbour[j_uprt] = FmmNeighbour()
            
            if i < k0:
                # Down-left neighbour (plus theta, minus phi)
                if th == no_theta-2:
                    j_dnlt = no_vertices - 1
                else:
                    if ph == 0:
                        j_dnlt = get_vertex_index(th+1, no_phi-1, no_vertices, no_theta, no_phi)
                    else:
                        j_dnlt = vertex_neighbour = get_vertex_index(th+1, ph-1, no_vertices, no_theta, no_phi)
                #vertex[i].neighbour.append(vertex_neighbour)
                vertex[i].neighbour[j_dnlt] = FmmNeighbour()
                
                # Down-right neighbour (plus theta, plus phi)
                if th == no_theta-2:
                    j_dnrt = no_vertices - 1
                else:
                    if ph == no_phi-1:
                        j_dnrt = get_vertex_index(th+1, 0, no_vertices, no_theta, no_phi)
                    else:
                        j_dnrt = get_vertex_index(th+1, ph+1, no_vertices, no_theta, no_phi)
                #vertex[i].neighbour.append(vertex_neighbour)
                vertex[i].neighbour[j_dnrt] = FmmNeighbour()
            
            # Find faces for 8-connectivity case
            if i <= no_phi:
                # North pole adjacent
                vertex[i].neighbour[j_up].face = [j_lt, j_rt]
                vertex[i].neighbour[j_up].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_up].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_rt].face = [j_up, j_dnrt]
                vertex[i].neighbour[j_rt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_rt].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face = [j_rt, j_dn]
                vertex[i].neighbour[j_dnrt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dn].face = [j_dnlt, j_dnrt]
                vertex[i].neighbour[j_dn].face_angle[j_dnlt] = cos_alpha
                vertex[i].neighbour[j_dn].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnlt].face = [j_dn, j_lt]
                vertex[i].neighbour[j_dnlt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dnlt].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_lt].face = [j_up, j_dnlt]
                vertex[i].neighbour[j_lt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_lt].face_angle[j_dnlt] = cos_alpha
                
            elif i >= k0:
                # South pole adjacent
                vertex[i].neighbour[j_up].face = [j_uplt, j_uprt]
                vertex[i].neighbour[j_up].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_up].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_uprt].face = [j_up, j_rt]
                vertex[i].neighbour[j_uprt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uprt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_rt].face = [j_uprt, j_dn]
                vertex[i].neighbour[j_rt].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_rt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dn].face = [j_lt, j_rt]
                vertex[i].neighbour[j_dn].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_dn].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_lt].face = [j_dn, j_uplt]
                vertex[i].neighbour[j_lt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_lt].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_uplt].face = [j_up, j_lt]
                vertex[i].neighbour[j_uplt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uplt].face_angle[j_lt] = cos_alpha
                
            else:
                # Not pole adjacent
                vertex[i].neighbour[j_up].face = [j_uplt, j_uprt]
                vertex[i].neighbour[j_up].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_up].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_uprt].face = [j_up, j_rt]
                vertex[i].neighbour[j_uprt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uprt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_rt].face = [j_uprt, j_dnrt]
                vertex[i].neighbour[j_rt].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_rt].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face = [j_rt, j_dn]
                vertex[i].neighbour[j_dnrt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dn].face = [j_dnlt, j_dnrt]
                vertex[i].neighbour[j_dn].face_angle[j_dnlt] = cos_alpha
                vertex[i].neighbour[j_dn].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnlt].face = [j_dn, j_lt]
                vertex[i].neighbour[j_dnlt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dnlt].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_lt].face = [j_uplt, j_dnlt]
                vertex[i].neighbour[j_lt].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_lt].face_angle[j_dnlt] = cos_alpha
                vertex[i].neighbour[j_uplt].face = [j_up, j_lt]
                vertex[i].neighbour[j_uplt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uplt].face_angle[j_lt] = cos_alpha
            
        else: 
            # Find faces for 4-connectivity case
            vertex[i].neighbour[j_up].face = [j_lt, j_rt]
            vertex[i].neighbour[j_up].face[j_lt] = cos_alpha
            vertex[i].neighbour[j_up].face[j_rt] = cos_alpha
            vertex[i].neighbour[j_dn].face = [j_lt, j_rt]
            vertex[i].neighbour[j_dn].face[j_lt] = cos_alpha
            vertex[i].neighbour[j_dn].face[j_rt] = cos_alpha
            vertex[i].neighbour[j_lt].face = [j_up, j_dn]
            vertex[i].neighbour[j_lt].face[j_up] = cos_alpha
            vertex[i].neighbour[j_lt].face[j_dn] = cos_alpha
            vertex[i].neighbour[j_rt].face = [j_up, j_dn]
            vertex[i].neighbour[j_rt].face[j_up] = cos_alpha
            vertex[i].neighbour[j_rt].face[j_dn] = cos_alpha

            
    #for j in vertex[i].neighbour:
    #    vertex[i].distance_to_neighbour[j] = -1.0

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
        
# Split pole facing angles
def split_polar_angle(vertex, i, w_pole, face_no, neighbour_no):
    cos_alpha = vertex[i].face_dot[face_no]
    j = vertex[vertex[i].neighbour[neighbour_no]].neighbour[neighbour_no]
    no_iterations = 0
    while cos_alpha < 0:
        w2 = vertex[j].carts - vertex[i].carts
        cos_alpha = np.dot(w_pole, w2)
        if cos_alpha >= 0:
            vertex[i].face.append([vertex[i].face[face_no][0], j])
            vertex[i].face.append([vertex[i].face[face_no][1], j])
            vertex[j].neighbour.append(i)
            break
        else:
            j = vertex[j].neighbour[neighbour_no]
            no_iterations += 1 
            if no_iterations == 360:
                print(['i = ', i, ', face ', face_no])
                raise NameError('No split could be found')
                        
# Split obtuse angles
def split_obtuse_angle(vertex, i, face_no, nb_orthog, nb_diag):
    w1 = vertex[vertex[i].neighbour[nb_orthog]].carts - vertex[i].carts
    cos_alpha = vertex[i].face_dot[face_no]
    j = vertex[vertex[i].neighbour[nb_diag]].neighbour[nb_orthog]
    no_iterations = 0
    while cos_alpha < 0:
        w2 = vertex[j].carts - vertex[i].carts
        cos_alpha = np.dot(w1, w2)
        if cos_alpha >= 0:
            vertex[i].face.append([vertex[i].face[face_no][0], j])
            vertex[i].face.append([vertex[i].face[face_no][1], j])
            vertex[j].neighbour.append(i)
            break
        else:
            j = vertex[j].neighbour[nb_orthog]
            no_iterations += 1 
            if no_iterations == 360:
                print(['i = ', i, ', face ', face_no])
                raise NameError('No split could be found')
                
    w3 = vertex[vertex[i].neighbour[nb_diag]].carts - vertex[i].carts
    
    if np.dot(w2, w3) < 0:
        print('invalid angle')
                    
# Split problematic triangles
def split_update_triangles(vertex, no_theta, no_phi, no_vertices):
    
    for i in range(no_vertices):
        
        th_index = get_theta_index(i, no_vertices, no_theta, no_phi)
        
        if th_index == 1: # North pole adjacent
            wp = vertex[0].carts - vertex[i].carts
            if vertex[i].face_valid[0] == False: # Split north-east triangle
                split_polar_angle(vertex, i, wp, 0, 3)
            if vertex[i].face_valid[1] == False: # Split east-south-east triangle
                split_obtuse_angle(vertex, i, 1, 3, 5)
            if vertex[i].face_valid[4] == False: # Split west-south-west triangle
                split_obtuse_angle(vertex, i, 4, 2, 4)
            if vertex[i].face_valid[5] == False: # Split north-west triangle
                split_polar_angle(vertex, i, wp, 5, 2)          
            no_faces = 6

        elif th_index == no_theta-2: # South pole adjacent
            wp = vertex[no_vertices-1].carts - vertex[i].carts
            if vertex[i].face_valid[1] == False: # Split east-north-east triangle
                split_obtuse_angle(vertex, i, 1, 3, 5)
            if vertex[i].face_valid[2] == False: # Split south-east triangle
                split_polar_angle(vertex, i, wp, 2, 3)       
            if vertex[i].face_valid[3] == False: # Split north-west triangle
                split_polar_angle(vertex, i, wp, 3, 2)
            if vertex[i].face_valid[4] == False: # Split west-north-west triangle
                split_obtuse_angle(vertex, i, 4, 2, 4)
            no_faces = 6        
            
                    
        else: # Not adjacent to a pole
            if vertex[i].face_valid[1] == False: # Split east-north-east triangle
                split_obtuse_angle(vertex, i, 1, 3, 5)
            if vertex[i].face_valid[2] == False: # Split east-south-east triangle
                split_obtuse_angle(vertex, i, 2, 3, 7)
            if vertex[i].face_valid[5] == False: # Split west-south-west triangle
                split_obtuse_angle(vertex, i, 5, 2, 6)
            if vertex[i].face_valid[6] == False: # Split west-north-west triangle
                split_obtuse_angle(vertex, i, 6, 2, 4)
            no_faces = 8
        
        m = 0        
        for j in range(no_faces):
            if vertex[i].face_valid[j] == False:
                del vertex[i].face[j-m]
                m += 1
                    
                    
               
            
            
            
# Check upward pointing, pole adjacent triangles
def find_obtuse_angles(vertex, no_theta, no_phi, no_vertices):
    no_obtuse = [0, 0, 0, 0, 0, 0, 0, 0]
    acute_north_pole = [[True]*no_phi]
    acute_north_pole_adjacent = [[True]*6] * no_phi
    acute_north = [[True]*8] * no_phi
    acute_equator = [[True]*8] * no_phi
    acute_south = [[True]*8] * no_phi
    acute_south_pole_adjacent = [[True]*6] * no_phi
    acute_south_pole = [[True]*no_phi]
    no_obtuse_north_adjacent = [0] * 6
    no_obtuse_north = [0] * 8
    no_obtuse_south = [0] * 8
    no_obtuse_south_adjacent = [0] * 6
    
    th_eq = (no_theta-1) / 2
    for i in range(no_vertices):
        if i == 0: # North pole
            for j in range(len(vertex[i].face_valid)):
                if vertex[i].face_valid[j] == False:
                    no_obtuse[0] += 1
                    no_obtuse[1] += 1
                    acute_north_pole[j] = False
        elif i == no_vertices-1: # South pole
            for j in range(len(vertex[i].face_valid)):
                if vertex[i].face_valid[j] == False:
                    no_obtuse[0] += 1
                    no_obtuse[7] += 1
                    acute_south_pole[j] = False
        else:
            th_index = get_theta_index(i, no_vertices, no_theta, no_phi)
            ph_index = get_phi_index(i, th_index, no_phi)
            if th_index == 1: # North pole adjacent
                for j in range(6):
                    if vertex[i].face_valid[j] == False:
                        no_obtuse[0] += 1
                        no_obtuse[2] += 1
                        acute_north_pole_adjacent[ph_index][j] = False
                        no_obtuse_north_adjacent[j] += 1
            elif th_index == no_theta-2: # South pole adjacent
                for j in range(6):
                    if vertex[i].face_valid[j] == False:
                        no_obtuse[0] += 1
                        no_obtuse[6] += 1
                        acute_south_pole_adjacent[ph_index][j] = False
                        no_obtuse_south_adjacent[j] += 1
            elif th_index == th_eq: # Equator
                for j in range(8):
                    if vertex[i].face_valid[j] == False:
                        no_obtuse[0] += 1
                        no_obtuse[4] += 1
                        acute_equator[ph_index][j] = False
            elif th_index < th_eq: # Northern hemisphere
                for j in range(8):
                    if vertex[i].face_valid[j] == False:
                        no_obtuse[0] += 1
                        no_obtuse[3] += 1
                        no_obtuse_north[j] += 1
                        if th_index == th_eq/2:
                            acute_north[ph_index][j] = False
            elif th_index > th_eq: # Southern hemisphere
                for j in range(8):
                    if vertex[i].face_valid[j] == False:
                        no_obtuse[0] += 1
                        no_obtuse[5] += 1
                        no_obtuse_south[j] += 1
                        if th_index == 3*th_eq/2:
                            acute_south[ph_index][j] = False
  
    return [no_obtuse, acute_north_pole, acute_north_pole_adjacent, acute_north,
            acute_equator, acute_south, acute_south_pole_adjacent, acute_south_pole,
            no_obtuse_north_adjacent, no_obtuse_north, no_obtuse_south, no_obtuse_south_adjacent]
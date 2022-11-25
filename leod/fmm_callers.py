# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:38:56 2022

@author: Callum Marples

Prepare and call the fast marching routine for a given precaluculated graph and
input start/end points.
"""

from numpy import array
from math import pi
from leod.fmm_polar_graph import find_vertex_index
from leod.fmm_fast_marching import fast_marching
from leod.fmm_fast_marching import endpoint_distance
from leod.triangulation import find_closest_vertex

def calculate_pair_distance(shape, vertex, start_point, end_point, order, graph_type="tri", grid=-1, is_radians=False):
    
    start_point_temp = [0.0, 0.0]
    end_point_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
    for i in range(2):
        start_point_temp[i] = start_point[i] * conv
        end_point_temp[i] = end_point[i] * conv
    
    # Get Cartesian coordinatess of start and end points.
    start_carts = shape.polar2cart(start_point_temp[0], start_point_temp[1])   
    end_carts = shape.polar2cart(end_point_temp[0], end_point_temp[1])
     
    # Find start and end vertices
    if graph_type == "polar":
        start_vertex, start_th_index, st_ph_index = find_vertex_index(grid.theta_list, grid.phi_list, start_point_temp[0], start_point_temp[1])
        end_vertex, end_th_index, end_ph_index = find_vertex_index(grid.theta_list, grid.phi_list, end_point_temp[0], end_point_temp[1])
    else:
        # Use triangulation
        start_vertex = find_closest_vertex(vertex, array(start_carts))
        end_vertex = find_closest_vertex(vertex, array(end_carts))
    
    
    end_vertex_carts = vertex[end_vertex].carts
    
    # If closest vertex and endpoint are the same, no interpolation needed.
    end_diff = end_vertex_carts - end_carts
    end_diff = end_diff[0]*end_diff[0] + end_diff[1]*end_diff[1] + end_diff[2]*end_diff[2]
    if end_diff > 1e-9:
        is_end_interpolate = True
    else:
        is_end_interpolate = False
    # Find neighbours of the endpoint, to use in the stopping criterion
    end_dict = {}
    end_dict[end_vertex] = False
    if is_end_interpolate == True:
        for j in vertex[end_vertex].neighbour.keys():
            end_dict[j] = False       
        
    # Call fast marching method
    fmm = fast_marching(vertex, start_vertex, start_carts, order, end_dict)
    
    if is_end_interpolate == True:
        d = endpoint_distance(vertex, fmm, end_carts, end_vertex, shape)
    else:
        d = fmm.distance[end_vertex]
    
    # Determine number of obtuse angles relevant to the endpoint distance
    temp_set = set()
    fmm.updater_set = set() # Set of vertices whose values influenced the endpoint
    for j in end_dict.keys():
        fmm.updater_set.add(j)
        temp_set.add(j)
    while len(temp_set) > 0:
        temp_set_new = set()
        for i in temp_set:
            for j in fmm.update[i]:
                if j != -1:
                    fmm.updater_set.add(j)
                    temp_set_new.add(j)
        temp_set = temp_set_new
    
    # Find obtuse vertices based on fmm.update
    fmm.no_obtuse_path = 0
    for i in fmm.updater_set:
        if len(fmm.update[i]) == 2 and fmm.update[i][1] == -1:
            fmm.no_obtuse_path += 1
    
    
    return d, fmm
        
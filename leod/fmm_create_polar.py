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
from .fmm_grid import FmmGrid

def generate_polar_grid(shape, no_theta, no_phi, is_Dijkstra=False, is_connect_8=False):
    
    # Create empty FmmGrid object
    grid = FmmGrid(is_Dijkstra)
    
    # Structured grid information
    no_pixels = (no_theta - 2)*no_phi + 2
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
    for i in range(no_pixels):
        
        # theta and phi indices
        th_index = get_theta_index(i, no_pixels, no_theta, no_phi)
        ph_index = get_phi_index(i, th_index, no_phi)
        polar_index.append([th_index, ph_index])
        
        # Cartesian coordinates
        carts = np.array( shape.polar2cart(theta_list[th_index], phi_list[ph_index]) )
        
        # Create new vertex
        grid.vertex.append(FmmVertex(i, carts))
        
    
    
    
    
    
    
    
    
    
    
    return grid

############################# Define subroutines #############################
    
# Get theta index from pixel index
def get_theta_index(pixel_index, no_pixels, no_theta, no_phi):
    if pixel_index > 0 and pixel_index < no_pixels-1:
        return math.ceil(float(pixel_index)/float(no_phi))
    elif pixel_index == 0:
        return 0
    elif pixel_index == no_pixels-1:
        return no_theta-1
    else:
        return no_theta
        
# Get phi index from pixel index and theta index
def get_phi_index(pixel_index, theta_index, no_phi):
    return ( pixel_index - 1 - no_phi*(theta_index-1) )
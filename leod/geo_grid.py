# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:30:04 2022

@author: Callum Marples
"""

import math
import numpy as np

from .binary_search import binary_search
from .ellipsoid_shape import EllipsoidShape
from .geo_pixel import GeoPixel

# This class contains the set of all GeoPixel objects for a given 
# EllipsoidShape. It also includes functions for initialising the set of 
# pixels, by calculating the neighbour-to-neighbour distances.
class GeoGrid:
    # Constructor
    def __init__(self, shape=EllipsoidShape(), no_theta=19, no_phi=36, is_refine=False):
        # Take inputs
        self.shape = shape
        self.no_theta = no_theta
        self.no_phi = no_phi
        self.is_refine = is_refine
        # Calculate number of pixels and required increments
        self.no_pixels = (self.no_theta - 2)*self.no_phi + 2
        self.delta_theta = math.pi / (self.no_theta - 1)
        self.delta_phi = 2.0*math.pi / self.no_phi
        # Compute lists of theta and phi values
        self.theta_list = [0.0] * self.no_theta
        self.phi_list = [0.0] * (self.no_phi + 1)
        for i in range(self.no_theta):
            self.theta_list[i] = i * self.delta_theta
        for i in range(self.no_phi+1):
            self.phi_list[i] = i * self.delta_phi
        
        if self.is_refine == True:
            self.is_initialised = False
            
        else:
            # Construct a list of pixel objects
            self.pixel = []
            for i in range(self.no_pixels):
                theta_index = self.get_theta_index(i)
                phi_index = self.get_phi_index(i, theta_index)
                carts = self.polars_to_cartesians( self.theta_list[theta_index],
                                                   self.phi_list[phi_index] )
                self.pixel.append(GeoPixel(i, theta_index, phi_index,
                                           carts, self.no_pixels) )
            # Set array sizes for pixel neighbour information
            self.find_neighbour_indices()
            self.is_initialised = False
            
    # Get pixel index from theta and phi indices
    def get_pixel_index(self, theta_index, phi_index):
        if theta_index > 0 and theta_index < self.no_theta-1:
            return 1 + phi_index + self.no_phi*(theta_index-1)
        elif theta_index == 0:
            return 0
        else: # theta_index = no_theta - 1
            return self.no_pixels - 1
        
    # Get theta index from pixel index
    def get_theta_index(self, pixel_index):
        if pixel_index > 0 and pixel_index < self.no_pixels-1:
            return math.ceil(float(pixel_index)/float(self.no_phi))
        elif pixel_index == 0:
            return 0
        elif pixel_index == self.no_pixels-1:
            return self.no_theta-1
        else:
            return self.no_theta
        
    # Get phi index from pixel index and theta index
    def get_phi_index(self, pixel_index, theta_index):
        return ( pixel_index - 1 - self.no_phi*(theta_index-1) )
    
    # Find theta index closest to given theta value, th
    def find_theta_index(self, th):
        return binary_search(self.theta_list, th)
    
    # Find phi index closest to given phi value, ph
    def find_phi_index(self, ph):
        index = binary_search(self.phi_list, ph)
        if index == self.no_phi:
            return 0
        else:
            return index
    
    # Calculate Cartesian coordinates of a given (theta, phi) point
    # on the ellipsoid, self.shape. Return the coordinates as a list
    def polars_to_cartesians(self, theta, phi):
        sin_theta = math.sin(theta)
        return [ self.shape.a_axis * sin_theta * math.cos(phi),
                 self.shape.b_axis * sin_theta * math.sin(phi),
                 self.shape.c_axis * math.cos(theta) ]
        
    # Find indices of neighbouring pixels
    # Each pixel has four neighbours. For non-polar pixels, neighbours are
    # always ordered as "Up, Down, Left, Right".
    def find_neighbour_indices(self):
        # Pole neighbours
        for i in range(self.no_phi):
            self.pixel[0].neighbour.append(1+i)
            self.pixel[0].neighbour_distance.append(-1.0)
            self.pixel[self.no_pixels-1].neighbour.append(self.no_pixels-2-i)
            self.pixel[self.no_pixels-1].neighbour_distance.append(-1.0)
        # Non-pole neighbours
        for i in range(1, self.no_theta-1):
            for j in range(self.no_phi):
                pixel_no = self.get_pixel_index(i, j)
                # Initialise distances
                self.pixel[pixel_no].neighbour_distance = [-1.0, -1.0, -1.0, -1.0]
                # Up neighbour (minus theta)
                if i == 1:
                    pixel_neighbour = 0
                else:
                    pixel_neighbour = self.get_pixel_index(i-1, j)
                self.pixel[pixel_no].neighbour.append(pixel_neighbour)
                # Down neighbour (plus theta)
                if i == self.no_theta-2:
                    pixel_neighbour = self.no_pixels - 1
                else:
                    pixel_neighbour = self.get_pixel_index(i+1, j)
                self.pixel[pixel_no].neighbour.append(pixel_neighbour)   
                # Left neighbour (minus phi)
                if j == 0:
                    pixel_neighbour = self.get_pixel_index(i, self.no_phi-1)
                else:
                    pixel_neighbour = self.get_pixel_index(i, j-1)
                self.pixel[pixel_no].neighbour.append(pixel_neighbour)
                # Right neighbour (plus phi)
                if j == self.no_phi-1:
                    pixel_neighbour = self.get_pixel_index(i, 0)
                else:
                    pixel_neighbour = self.get_pixel_index(i, j+1)
                self.pixel[pixel_no].neighbour.append(pixel_neighbour) 
                
    # Get Euclidean distance between pixel i and neighbour k.
    # If this distance has not been calculated yet, then do so.
    def get_distance(self, pix_i, k):
        # pix_i : first pixel
        # k : neighbour number (if non-pole, must be 1, 2, 3, or 4)
        # j : index of second pixel (calculated from i and k)
        if pix_i.neighbour_distance[k] == -1.0:
            # Calculate Euclidean distance between pixel i and its neighbour j
            j = pix_i.neighbour[k]
            pix_i.neighbour_distance[k] = math.sqrt( np.sum((pix_i.carts - self.pixel[j].carts)**2.0) )
            # Assign this distance to the relevant neighbour of pixel j.
            # Neighbours are always ordered as "0, 1, 2, 3" = "Up, Down, Left, Right".
            # e.g. if j is the "down" neighbour of i, then i must be the "up" neighbour of j.
            if pix_i.theta_index == 0:
                self.pixel[j].neighbour_distance[0] = pix_i.neighbour_distance[k]
            elif k == 0 or pix_i.theta_index == self.no_theta-1:
                self.pixel[j].neighbour_distance[1] = pix_i.neighbour_distance[k]
            elif k == 1:
                self.pixel[j].neighbour_distance[0] = pix_i.neighbour_distance[k]
            elif k == 2:
                self.pixel[j].neighbour_distance[3] = pix_i.neighbour_distance[k]
            else:
                self.pixel[j].neighbour_distance[2] = pix_i.neighbour_distance[k]
        return pix_i.neighbour_distance[k]
    
    # Initialise a grid to be used in a source refinement fast marching calculation
    def initialise_refined_grid(self, no_theta_border, no_phi_border, centre_theta, centre_phi):
        print('rfnd')
        
    
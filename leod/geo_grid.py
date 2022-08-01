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
                self.pixel.append(GeoPixel(i, theta_index, phi_index, self.no_pixels,
                                           carts, False) )
            # Set array sizes for pixel neighbour information
            for i in range(self.no_pixels):    
                th = self.get_theta_index(i)
                ph = self.get_phi_index(i, th)
                self.find_neighbour_indices(i, th, ph)
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
    
    # Find index of pixel containing given theta and phi
    def find_pixel_index(self, th, ph):
        th_index = self.find_theta_index(th)
        ph_index = self.find_phi_index(ph)
        return self.get_pixel_index(th_index, ph_index)
        
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
    def find_neighbour_indices(self, pixel_no, th, ph):
        # Pole neighbours
        if pixel_no == 0:                  # North pole
            for k in range(self.no_phi):
                self.pixel[0].neighbour.append(1+k)
                self.pixel[0].neighbour_distance.append(-1.0)
        elif pixel_no == self.no_pixels-1: # South pole
            for k in range(self.no_phi):
                self.pixel[self.no_pixels-1].neighbour.append(self.no_pixels-2-k)
                self.pixel[self.no_pixels-1].neighbour_distance.append(-1.0)
        else:                              # Not a pole
            # Initialise distances
            self.pixel[pixel_no].neighbour_distance = [-1.0, -1.0, -1.0, -1.0]
            # Up neighbour (minus theta)
            if th == 1:
                pixel_neighbour = 0
            else:
                pixel_neighbour = self.get_pixel_index(th-1, ph)
            self.pixel[pixel_no].neighbour.append(pixel_neighbour)
            # Down neighbour (plus theta)
            if th == self.no_theta-2:
                pixel_neighbour = self.no_pixels - 1
            else:
                pixel_neighbour = self.get_pixel_index(th+1, ph)
            self.pixel[pixel_no].neighbour.append(pixel_neighbour)   
            # Left neighbour (minus phi)
            if ph == 0:
                pixel_neighbour = self.get_pixel_index(th, self.no_phi-1)
            else:
                pixel_neighbour = self.get_pixel_index(th, ph-1)
            self.pixel[pixel_no].neighbour.append(pixel_neighbour)
            # Right neighbour (plus phi)
            if ph == self.no_phi-1:
                pixel_neighbour = self.get_pixel_index(th, 0)
            else:
                pixel_neighbour = self.get_pixel_index(th, ph+1)
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
    def initialise_refined_grid(self, no_theta_border, no_phi_border, centre_theta, centre_phi, pix_refine, main2rfnd):
        
        
        
        # Check if north pole is within the border
        if centre_theta == no_theta_border or (centre_theta != 0 and centre_theta < no_theta_border):
            self.border_pixels.append(0)               # North pole on border
            self.pixel.append(GeoPixel(0, 0, 0, self.no_pixels, np.array([0,0,self.shape.c_axis]), False))
            self.find_neighbour_indices(0, 0, 0)
        else:
            self.pixel.append(GeoPixel(0, 0, 0, self.no_pixels, [0,0,0], True))
        
        # Check all pixels for inclusion in the refined grid
        # Only include pixel if it is within a main pixel in the refinement range
        for i in range(1, self.no_pixels-1):
            th = self.get_theta_index(i)
            ph = self.get_phi_index(i, th)
            pixel_valid = True
            if centre_theta == 0:                 # Centre pixel is north pole
                if th > no_theta_border:
                    pixel_valid == False
                elif th == no_theta_border:
                    self.border_pixels.append(i)
            elif centre_theta == self.no_theta-1: # Centre pixel is south pole
                if th < self.no_theta - 1 - no_theta_border:
                    pixel_valid == False
                elif th == self.no_theta - 1 - no_theta_border:
                    self.border_pixels.append(i)
            else:                                 # Centre pixel not a pole
                th_diff = abs(th - centre_theta)
                ph_diff = ph - centre_phi
                # Handle phi periodicity
                if ph_diff > self.no_phi/2:
                    ph_diff -= self.no_phi
                elif ph_diff < -self.no_phi/2:
                    ph_diff += self.no_phi
                ph_diff = abs(ph_diff)
                # Check for inclusion
                if th_diff > no_theta_border or abs(ph_diff) > no_phi_border:
                    pixel_valid = False
                # Check for border
                if ((pixel_valid == True) and 
                   (th_diff == no_theta_border or abs(ph_diff) == no_phi_border)):
                   self.border_pixels.append(i)
            if pixel_valid == True:
                # Initialise full GeoPixel
                carts = self.polars_to_cartesians( self.theta_list[th],
                                                   self.phi_list[ph] )
                self.pixel.append(GeoPixel(i, th, ph, self.no_pixels,
                                           carts, False) )
                self.find_neighbour_indices(i, th, ph)
            else:
                # Initialise GeoPixel, containing only the index information
                self.pixel.append(GeoPixel(i, th, ph, self.no_pixels, [0,0,0], True))
                
        # Check if south pole is within the border        
        if ((self.no_theta - 1 - centre_theta == no_theta_border) or 
             (centre_theta != self.no_theta - 1 and self.no_theta - 1 - centre_theta < no_theta_border)):
            self.border_pixels.append(self.no_pixels-1) # South pole on border
            self.pixel.append(GeoPixel(self.no_pixels-1, self.no_theta-1, 0, self.no_pixels,
                                       np.array([0,0,-self.shape.c_axis]), False))
            self.find_neighbour_indices(self.no_pixels-1, self.no_theta-1, 0)
        else:
            self.pixel.append(GeoPixel(self.no_pixels-1, self.no_theta-1, 0,
                                       self.no_pixels, [0,0,0], True))
        
        
               
                










        
    
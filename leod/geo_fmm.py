# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:16:23 2022

@author: Callum Marples
"""

# In this class, the distances from a starting point to all other pixels in a
# GeoGrid object (or to a subset of them) are approximated using the Fast 
# Marching Method.

import math
import numpy as np
import heapq
from .geo_pixel import GeoPixel

# Calculate the Euclidean distance between the representative points
# of pixels p and q.
def euclidean_distance(p, q):
    return math.sqrt( np.sum((p.carts - q.carts)**2.0) )


class GeoFMM:
    # Constructor
    def __init__(self, grid, theta, phi):
        self.theta = theta
        self.phi = phi
        self.grid = grid
        self.initialise(theta, phi)
        
    # Initialise fast marching method for a given start point
    def initialise(self, theta, phi):
        # Set Start pixel
        self.start_pixel = self.initialise_pixel(theta, phi) 
        # Initialise distance from the start and set all pixels to not alive
        self.geo_distances = [math.inf] * self.grid.no_pixels
        self.alive = [False] * self.grid.no_pixels

        
        
    # Initialise neighbours and distances for start/end pixel
    def initialise_pixel(self, theta, phi):
        theta_index = self.grid.find_theta_index(theta)
        phi_index =self. grid.find_phi_index(phi)
        pixel_index = self.grid.get_pixel_index(theta_index, phi_index)
        carts = self.grid.polars_to_cartesians(theta, phi)
        pix = GeoPixel(pixel_index, theta_index, phi_index,
                       carts, self.grid.no_pixels)
        pix.neighbour = []
        pix.neighbour_distance = []
        if pix.is_north == True:   # North pole
            pix.phi_index = self.grid.no_phi - 1
            for i in range(self.grid.no_phi):
                neigh = 1 + i
                pix.neighbour.append(neigh)
                pix.neighbour_distance.append( euclidean_distance(pix, self.grid.pixel[neigh]) )
        elif pix.is_south == True: # South pole
            pix.phi_index = self.grid.no_phi - 1
            for i in range(self.grid.no_phi):
                neigh = self.grid.no_pixels - 2 - i
                pix.neighbour.append(neigh)
                pix.neighbour_distance.append( euclidean_distance(pix, self.grid.pixel[neigh]) )
        else:                      # Not a pole
            i = pix.theta_index
            j = pix.phi_index
            # Up neighbour (minus theta)
            if i == 1:
                pixel_neighbour = 0
            else:
                pixel_neighbour = self.grid.get_pixel_index(i-1, j)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, self.grid.pixel[pixel_neighbour]) )
            # Down neighbour (plus theta)
            if i == self.grid.no_theta-2:
                pixel_neighbour = self.grid.no_pixels - 1
            else:
                pixel_neighbour = self.grid.get_pixel_index(i+1, j)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, self.grid.pixel[pixel_neighbour]) )
            # Left neighbour (minus phi)
            if j == 0:
                pixel_neighbour = self.grid.get_pixel_index(i, self.no_phi-1)
            else:
                pixel_neighbour = self.grid.get_pixel_index(i, j-1)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, self.grid.pixel[pixel_neighbour]) )
            # Right neighbour (plus phi)
            if j == self.grid.no_phi-1:
                pixel_neighbour = self.grid.get_pixel_index(i, 0)
            else:
                pixel_neighbour = self.grid.get_pixel_index(i, j+1)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, self.grid.pixel[pixel_neighbour]) )
        return pix
    
    # Calculate geodesic distances using the fast marching method, or Dijkstra's algorithm.
    # If a single endpoint is desired, input positive values of theta_end and
    # phi_end. If either is negative, then distance to all pixels is found.
    def calculate_geodesics(self, order, theta_end=-1.0, phi_end=-1.0):
        if theta_end >= 0.0 and phi_end >= 0.0:
            # Compute shortest distance from start point to end point.
            # Set end pixel.
            self.end_pixel = self.initialise_pixel(theta_end, phi_end)
            # Perfrom wavefront propagation loop.
            
            return -1.0
        else:
            # Compute shortest distance from start point to all other pixels.
            # Perfrom wavefront propagation loop.
            
            return -1.0
                
                
    # Determine the required shortest distances by wavefront propagation.
    # This is the main loop for Dijkstra's algorithm and the fast marching method.            
    def perform_loop(self, order, end_flag, refine_flag):
        if end_flag == True:
            self.end_flag = True
        else:
            self.end_flag == False
        queue = [(0.0, self.start_pixel.pixel_index)]
        no_alive = 0 # Number of pixels with a set distance
        # Perform wavefront propagation
        while no_alive != self.grid.no_pixels:
            # Obtain index of pixel with smallest distance
            trial_distance, trial = heapq.heappop(queue)
            # If obtained trial index already alive, skip iteration
            if self.alive[trial] == True:
                continue
            # Set trial pixel to alive
            self.alive[trial] = True
            no_alive += 1
            # Check end point
            if end_flag == True:
                if trial == self.end_pixel.pixel_index:
                    # end_distance = ...
                    break
            # Visit neighbours of trial pixel
            for i in range(len(self.grid.pixel[trial].neighbour)):
                visit = self.grid.pixel[trial].neighbour[i]
                if self.alive[visit] == True:
                    continue
                #proposed_distance = ...
        
        
        
    # Determine distance from a trial to visited pixel
    def fmm_distance(self, trial, visit, neighbour_no, order):
        if order > 0: # Fast marching method
            proposed_distance = -1.0
        else:         # Dijkstra's algorithm
            proposed_distance = self.point_to_point_distance(trial, neighbour_no)
        return proposed_distance
    
    # Calculate distance between representative points of trial and neighbour
    def point_to_point_distance(self, trial, neighbour_no):
        if trial == self.start_pixel.pixel_index:
            proposed_distance = self.start_pixel.neighbour_distance[neighbour_no]
        else:
            if self.end_flag == True:
                if self.grid.pixel[trial].neighbour[neighbour_no] == self.end_pixel.pixel_index:
                    proposed_distance = self.neighbour_distance(trial, neighbour_no)
                else:
                    proposed_distance = self.grid.get_distance(self.grid.pixel[trial], neighbour_no)
            else:
                proposed_distance = self.grid.get_distance(self.grid.pixel[trial], neighbour_no)
        proposed_distance += self.geo_distances[trial]
        return proposed_distance
        
        
        
        
        
        return 0
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
                pixel_neighbour = self.grid.get_pixel_index(i, self.grid.no_phi-1)
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
    def calculate_geodesics(self, order, theta_end=-1.0, phi_end=-1.0,
                            is_refine=False, refine_range=1, refine_theta=3, refine_phi=3):
        if is_refine == True:
            self.refine_range = refine_range
            self.refine_theta = refine_theta
            self.refine_phi = refine_phi
            self.initialise_refined_grid()
            return -1.0
        else:
            if theta_end >= 0.0 and phi_end >= 0.0:
                # Compute shortest distance from start point to end point.
                # Set end pixel.
                self.end_pixel = self.initialise_pixel(theta_end, phi_end)
                # Perfrom wavefront propagation loop.
                self.perform_loop(order, end_flag=True, refine_flag=False)
                return self.end_distance
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
        self.geo_distances[self.start_pixel.pixel_index] = 0.0
        queue = [(0.0, self.start_pixel.pixel_index)]
        no_alive = 0 # Number of pixels with a set distance
        end_prev = -1
        end_neighbour = -1
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
            if self.end_flag == True:
                if trial == self.end_pixel.pixel_index:
                    self.end_distance = self.fmm_distance(end_prev, self.end_pixel.pixel_index, end_neighbour, order)
                    break
                
            # Visit neighbours of trial pixel
            for i in range(len(self.grid.pixel[trial].neighbour)):
                visit = self.grid.pixel[trial].neighbour[i]
                if self.alive[visit] == True:
                    continue
                proposed_distance = self.fmm_distance(trial, visit, i, order)
                if proposed_distance < self.geo_distances[visit]:
                    # Add visited pixel to the queue and update value in geo_distances
                    self.geo_distances[visit] = proposed_distance
                    heapq.heappush(queue, (proposed_distance, visit))
                    # If appropriate, record pixel used to derive this distance
                    if self.end_flag == True:
                        if visit == self.end_pixel.pixel_index:
                            end_prev = trial
                            end_neighbour = i
                        
        
        
        
    # Determine distance from a trial to visited pixel
    def fmm_distance(self, trial, visit, neighbour_no, order):
        if order > 0: # Fast marching method
            if visit == 0 or visit == self.grid.no_pixels-1: # If visited pixel is a pole
                proposed_distance = -1.0
            else:
                proposed_distance = self.fast_marching_method(visit, order)
            # If the fast marching method is not applicable, get a distance of -1
            # and Dijkstra's algorithm (point-to-point distance) is used instead.
            if proposed_distance == -1.0:
                proposed_distance = self.point_to_point_distance(trial, neighbour_no)
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
    
    # Determine neighbour-to-neighbour distance (including start/end pixel cases)
    def neighbour_distance(self, visit, neighbour_no):
        neighbour_dist = 0.0
        if self.grid.pixel[visit].neighbour[neighbour_no] == self.start_pixel.pixel_index:
            # Get neighbour number of visit pixel relative to start pixel
            if self.grid.pixel[visit].neighbour[neighbour_no] == 0: # Neighbour is north pole
                neighbour_dist = self.start_pixel.neighbour_distance[visit-1]
            elif self.grid.pixel[visit].neighbour[neighbour_no] == self.grid.no_pixels-1: # Neighbour is south pole
                neighbour_dist = self.start_pixel.neighbour_distance[self.grid.no_pixels-2-visit]
            else: # Neighbour not a pole
                if neighbour_no == 0:
                    neighbour_dist = self.start_pixel.neighbour_distance[1]
                elif neighbour_no == 1:
                    neighbour_dist = self.start_pixel.neighbour_distance[0]
                elif neighbour_no == 2:
                    neighbour_dist = self.start_pixel.neighbour_distance[3]
                elif neighbour_no == 3:
                    neighbour_dist = self.start_pixel.neighbour_distance[2]
        elif self.end_flag == True:
            if visit == self.end_pixel.pixel_index:
                neighbour_dist = self.end_pixel.neighbour_distance[neighbour_no]
            else:
                neighbour_dist = self.grid.get_distance(self.grid.pixel[visit], neighbour_no)
        else:
            neighbour_dist = self.grid.get_distance(self.grid.pixel[visit], neighbour_no)
        return neighbour_dist
    
    
    
    # Calculate distances using the fast marching method
    def fast_marching_method(self, visit, order):
        # Check for existence of 'alive' neighbours
        if ( ( self.alive[self.grid.pixel[visit].neighbour[0]] == True or self.alive[self.grid.pixel[visit].neighbour[1]] == True  ) and
             ( self.alive[self.grid.pixel[visit].neighbour[2]] == True or self.alive[self.grid.pixel[visit].neighbour[3]] == True  ) ):
            # i.e. If up or down neighbour alive AND left or right neighbour alive, use fast marching method
            
            # Determine whether Up or Down neighbour has smallest distance from start - use for theta distance.
            pix_theta = self.grid.pixel[visit].neighbour[0]
            neighbour_theta = 0     # Up
            dist_theta = self.geo_distances[pix_theta]
            if self.geo_distances[self.grid.pixel[visit].neighbour[1]] < dist_theta:
                pix_theta = self.grid.pixel[visit].neighbour[1]
                dist_theta = self.geo_distances[pix_theta]
                neighbour_theta = 1 # Down
            
            # Determine whether Left or Right neighbour has smallest distance from start - use for phi distance. 
            pix_phi = self.grid.pixel[visit].neighbour[2]
            neighbour_phi = 2     # Left
            dist_phi = self.geo_distances[pix_phi]
            if self.geo_distances[self.grid.pixel[visit].neighbour[3]] < dist_phi:
                pix_phi = self.grid.pixel[visit].neighbour[3]
                dist_phi = self.geo_distances[pix_phi]
                neighbour_phi = 3 # Right
            
            # Get distances from visited pixel to these neighbours
            delta_theta = self.neighbour_distance(visit, neighbour_theta)
            delta_phi = self.neighbour_distance(visit, neighbour_phi)
            
            # Use the fast marching method to calculate distance from start to visited pixel
            if order == 2:
                # 2nd Order FMM
                
                # Check for second order theta neighbours
                pole_flag = True
                pix_theta_2 = -1
                dist_theta_2 = dist_theta + 1
                if pix_theta != 0 and pix_theta != self.grid.no_pixels-1: # i.e. pix_theta is not a pole
                    pix_theta_2 = self.grid.pixel[pix_theta].neighbour[neighbour_theta]
                    dist_theta_2 = self.geo_distances[pix_theta_2]
                    pole_flag = False
                
                # theta terms of quadratic
                if pole_flag == False and self.alive[pix_theta_2] == True and dist_theta_2 <= dist_theta:
                    delta_theta_2 = self.neighbour_distance(pix_theta, neighbour_theta)
                    terms_theta = self.quadratic_terms(dist_theta, dist_theta_2, delta_theta, delta_theta_2, 2)
                else:
                    terms_theta = self.quadratic_terms(dist_theta, 0.0, delta_theta, 0.0, 1)
                    
                # Check for second order phi neighbours
                pix_phi_2 = self.grid.pixel[pix_phi].neighbour[neighbour_phi]
                dist_phi_2 = self.geo_distances[pix_phi_2]
                
                # Phi terms of quadratic
                if self.alive[pix_phi_2] and dist_phi_2 <= dist_phi:
                    delta_phi_2 = self.neighbour_distance(pix_phi, neighbour_phi)
                    terms_phi = self.quadratic_terms(dist_phi, dist_phi_2, delta_phi, delta_phi_2, 2)
                else:
                    terms_phi = self.quadratic_terms(dist_phi, 0.0, delta_phi, 0.0, 1)
                
                # Combine terms and subtract one from constant term
                terms = terms_theta + terms_phi
                terms[2] -= 1.0
                # Solve quadratic for distance from start point
                proposed_distance = self.solve_quadratic(terms)
                
            else:
                # 1st Order FMM
             
                # Square reciprocals of the visit-to-neighbour distances
                r2_theta = 1.0 / (delta_theta*delta_theta)
                r2_phi = 1.0 / (delta_phi*delta_phi)
             
                # Set up quadratic
                terms = [0] * 3 # [alpha, beta, gamma] with quadratic alpha*x^2 + beta*x + gamma = 0
                terms[0] = r2_theta + r2_phi
                terms[1] = -2.0*(dist_theta*r2_theta + dist_phi*r2_phi)
                terms[2] = (dist_theta/delta_theta)**2 + (dist_phi/delta_phi)**2 - 1.0
                proposed_distance = self.solve_quadratic(terms)
        else: # If no alive neighbour exists, use point-to-point distance
            proposed_distance = -1.0
        return proposed_distance
                
                    
    
    # Calculate quadratic terms for the fast marching method.
    # This only handles one of two directions (theta or phi) and is used
    # in the 2nd order FMM when only one double neighbour exists.
    def quadratic_terms(self, dist1, dist2, delta1, delta2, order):
        terms = np.zeros(3) # [alpha, beta, gamma] with quadratic alpha*x^2 + beta*x + gamma = 0
        if order == 2:
            terms[0] = 9.0 / (4.0*delta1*delta1)
            terms[1] = -3.0*(dist1 - dist2)/(2.0*delta1*delta2) - 2.0*terms[0]*dist1
            terms[2] = -3.0*dist1*(dist2 - dist1)/(2.0*delta1*delta2) + terms[0]*dist1*dist1 + (dist1*dist1 + dist2*dist2 - 2.0*dist1*dist2)/(4.0*delta2*delta2)
        elif order == 1:
            terms[0] = 1.0/(delta1*delta1)
            terms[1] = -2.0*terms[0]*dist1
            terms[2] = terms[0]*dist1*dist1
        return terms
    
    # Solve quadratic for fast marching method, using terms obtained from the
    # quadratic_terms function.
    def solve_quadratic(self, terms):
        discr = terms[1]*terms[1] - 4.0*terms[0]*terms[2]
        if discr > 0:
            proposed_distance = (-terms[1] + math.sqrt(discr)) / (2.0*terms[0])
        else:
            proposed_distance = -1.0 # Use point-to-point distance via Dijkstra's algorithm instead.
        return proposed_distance
    
    # Initialise a refined grid for use in a source refined fast marching method.
    def initialise_refined_grid(self):
        # pix_refine is a list containing all main grid pixels within the 
        # refinement range of the start pixel.
        self.pix_refine = []
        if self.start_pixel.pixel_index == 0:                       # North pole
            self.pix_refine.append(0)
            for i in range(self.refine_range):
                for j in range(self.grid.n_phi):
                    self.pix_refine.append(self.grid.get_pixel_index(i+1, j))
        elif self.start_pixel.pixel_index == self.grid.no_pixels-1: # South pole
            self.pix_refine.append(self.grid.no_pixels-1)
            for i in range(self.refine_range):
                for j in range(self.grid.n_phi):
                    self.pix_refine.append(self.grid.get_pixel_index(self.grid.no_theta - 2 - i, j))
        else:                                                       # Non-pole
             for i in range(self.grid.no_pixels):
                 th = self.grid.pixel[i].theta_index
                 ph = self.grid.pixel[i].phi_index
                 if ( (abs(th - self.start_pixel.theta_index) <= self.refine_range) and
                    (  abs(ph - self.start_pixel.phi_index) <= self.refine_range or
                       abs(ph - self.start_pixel.phi_index) >= self.grid.no_phi - self.refine_range or
                       self.grid.pixel[i].is_north == True or self.grid.pixel[i].is_south == True) ):
                     self.pix_refine.append(i)
        # Create refined grid
        no_theta_rfnd = self.grid.no_theta * self.refine_theta - 2
        no_phi_rfnd = self.grid.no_phi * self.refine_phi
        no_theta_border = int(self.refine_range*self.refine_theta + (self.refine_theta-1)/2)
        no_phi_border = int(self.refine_range*self.refine_phi + (self.refine_phi-1)/2)
         
         
                
                
                

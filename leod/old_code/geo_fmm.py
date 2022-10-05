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
import copy
from .geo_pixel import GeoPixel
from .geo_grid import GeoGrid

# Calculate the Euclidean distance between the representative points
# of pixels p and q.
def euclidean_distance(p, q):
    return math.sqrt( np.sum((p.carts - q.carts)**2.0) )


class GeoFMM:
    # Constructor
    def __init__(self, grid, theta, phi, is_flat=False):
        self.theta = theta
        self.phi = phi
        self.grid_main = grid
        self.is_flat = is_flat
        self.initialise(theta, phi, grid)
        
    # Initialise fast marching method for a given start point
    def initialise(self, theta, phi, grid):
        # Set Start pixel
        self.start_pixel = self.initialise_pixel(theta, phi, grid, is_end=False) 
        self.is_refined_grid = False
        # Initialise distance and alive lists as empty
        self.geo_distances = {}
        self.alive = {}
        # Set inverted list of neighbours
        if grid.is_neighbour8 == True:
            self.inverted_neighbour = [1, 0, 3, 2, 7, 6, 5, 4]
        else:
            self.inverted_neighbour = [1, 0, 3, 2]

        
        
    # Initialise neighbours and distances for start/end pixel
    def initialise_pixel(self, theta, phi, grid, is_end=False):
        theta_index = grid.find_theta_index(theta)
        phi_index = grid.find_phi_index(phi)
        pixel_index = grid.get_pixel_index(theta_index, phi_index)
        if self.is_flat == False:
            carts = grid.polars_to_cartesians(theta, phi)
        else:
            carts = [phi, theta]
        pix = GeoPixel(pixel_index, theta_index, phi_index,  grid.no_pixels,
                       carts, is_endpoint=is_end, is_refine=False)
        pix.neighbour = []
        pix.neighbour_distance = []
        if pix.is_north == True:   # North pole
            pix.phi_index = grid.no_phi - 1
            for i in range(grid.no_phi):
                neigh = 1 + i
                pix.neighbour.append(neigh)
                pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[neigh]) )
        elif pix.is_south == True: # South pole
            pix.phi_index = grid.no_phi - 1
            for i in range(grid.no_phi):
                neigh = grid.no_pixels - 2 - i
                pix.neighbour.append(neigh)
                pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[neigh]) )
        else:                      # Not a pole
            i = pix.theta_index
            j = pix.phi_index
            # Up neighbour (minus theta)
            if i == 1:
                pixel_neighbour = 0
            else:
                pixel_neighbour = grid.get_pixel_index(i-1, j)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
            # Down neighbour (plus theta)
            if i == grid.no_theta-2:
                pixel_neighbour = grid.no_pixels - 1
            else:
                pixel_neighbour = grid.get_pixel_index(i+1, j)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
            # Left neighbour (minus phi)
            if j == 0:
                pixel_neighbour = grid.get_pixel_index(i, grid.no_phi-1)
            else:
                pixel_neighbour = grid.get_pixel_index(i, j-1)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
            # Right neighbour (plus phi)
            if j == grid.no_phi-1:
                pixel_neighbour = grid.get_pixel_index(i, 0)
            else:
                pixel_neighbour = grid.get_pixel_index(i, j+1)
            pix.neighbour.append(pixel_neighbour)
            pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
            if grid.is_neighbour8 == True:
                # Up-left neighbour (minus theta, minus phi)
                if i == 1:
                    pixel_neighbour = 0
                else:
                    if j == 0:
                        pixel_neighbour = grid.get_pixel_index(i-1, grid.no_phi-1)
                    else:
                        pixel_neighbour = grid.get_pixel_index(i-1, j-1)
                pix.neighbour.append(pixel_neighbour)
                pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
                # Up-right neighbour (minus theta, plus phi)
                if i == 1:
                    pixel_neighbour = 0
                else:
                    if j == grid.no_phi-1:
                        pixel_neighbour = grid.get_pixel_index(i-1, 0)
                    else:
                        pixel_neighbour = grid.get_pixel_index(i-1, j+1)
                pix.neighbour.append(pixel_neighbour)
                pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
                # Down-left neighbour (plus theta, minus phi)
                if i == grid.no_theta-2:
                    pixel_neighbour = grid.no_pixels - 1
                else:
                    if j == 0:
                        pixel_neighbour = grid.get_pixel_index(i+1, grid.no_phi-1)
                    else:
                        pixel_neighbour = pixel_neighbour = grid.get_pixel_index(i+1, j-1)
                pix.neighbour.append(pixel_neighbour)
                pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
                # Up-right neighbour (minus theta, plus phi)
                if i == grid.no_theta-2:
                    pixel_neighbour = grid.no_pixels - 1
                else:
                    if j == grid.no_phi-1:
                        pixel_neighbour = grid.get_pixel_index(i+1, 0)
                    else:
                        pixel_neighbour = grid.get_pixel_index(i+1, j+1)
                pix.neighbour.append(pixel_neighbour)
                pix.neighbour_distance.append( euclidean_distance(pix, grid.pixel[pixel_neighbour]) )
        return pix
    
    
    # Set list of end pixels and a dictionary to determine whether endpoints exist
    # and which ones to use
    def set_end_pixels(self, grid, is_refine=False):
        self.end_pixel = []
        self.end_pixel_dict = {}
        self.end_distance = []
        self.no_endpoints_alive = 0
        for i in range(len(self.theta_end)):
            end_index = grid.find_pixel_index(self.theta_end[i], self.phi_end[i])
            if is_refine == False or end_index in self.main2rfnd:
                self.end_pixel.append(self.initialise_pixel(self.theta_end[i], self.phi_end[i], grid, is_end=True))
            else:
                self.end_pixel.append(-1)
            if end_index in self.end_pixel_dict:
                self.end_pixel_dict[end_index].append(i)
            else:
                self.end_pixel_dict[end_index] = []
                self.end_pixel_dict[end_index].append(i)
            self.end_distance.append(-1.0)
    
    # Transfer endpoint pixels from the refined grid, if that end distance has not yet been found
    def transfer_end_pixels(self, grid):
        self.end_refined = []
        self.end_pixel_dict = {}        
        for i in range(len(self.theta_end)):
            self.end_refined.append(False)
            end_index_main = grid.find_pixel_index(self.theta_end[i], self.phi_end[i])
            if self.end_pixel[i] == -1:
                self.end_pixel[i] = self.initialise_pixel(self.theta_end[i], self.phi_end[i], grid, is_end=True)
            else: # Refined grid GeoPixel exists
                end_index_rfnd = self.end_pixel[i].pixel_index
                if self.alive[end_index_rfnd] == False:
                    self.end_pixel[i] = self.initialise_pixel(self.theta_end[i], self.phi_end[i], grid, is_end=True)
                else: # This endpoint has already been found
                    self.end_refined[i] = True
            if end_index_main in self.end_pixel_dict:
                self.end_pixel_dict[end_index_main].append(i)
            else:
                self.end_pixel_dict[end_index_main] = []
                self.end_pixel_dict[end_index_main].append(i)    
                
    # Calculate geodesic distances using the fast marching method, or Dijkstra's algorithm.
    # If a single endpoint is desired, input positive values of theta_end and
    # phi_end. If either is negative, then distance to all pixels is found.
    def calculate_geodesics(self, order, theta_end=[], phi_end=[],
                            is_refine=False, refine_range=10, refine_theta=5, refine_phi=5):
        
        self.theta_end = theta_end
        self.phi_end = phi_end
        
        if is_refine == True:
            
            self.refine_range = refine_range
            self.refine_theta = refine_theta
            self.refine_phi = refine_phi
            self.initialise_refined_grid(self.grid_main)
            
            if len(theta_end) > 0 and len(phi_end) > 0:
                
                assert len(theta_end) == len(phi_end)
                self.set_end_pixels(self.grid_rfnd, is_refine=True)
                
                self.is_refined_grid = True
                self.perform_loop(order, self.grid_rfnd, self.start_pixel_rfnd, end_flag=False, refine_flag=True, main_flag=False)
                if -1.0 in self.end_distance:
                    self.transfer_end_pixels(self.grid_main)
                    self.transfer_grid()
                    self.is_refined_grid = False
                    self.perform_loop(order, self.grid_main, self.start_pixel, end_flag=True, refine_flag=False, main_flag=True)
                return self.end_distance
                
            else:
                return -1.0
        else:
            
            for i in range(self.grid_main.no_pixels):
                self.alive[i] = False
                self.geo_distances[i] = math.inf
            
            if len(theta_end) > 0 and len(phi_end) > 0:
                # Compute shortest distance from start point to end points.
                assert len(theta_end) == len(phi_end)
                # Set end pixels.
                self.set_end_pixels(self.grid_main)
                
                # Perfrom wavefront propagation loop.
                self.perform_loop(order, self.grid_main, self.start_pixel, end_flag=True)
                return self.end_distance
            else:
                # Compute shortest distance from start point to all other pixels.
                # Perfrom wavefront propagation loop.
                
                return -1.0
                
                
    # Determine the required shortest distances by wavefront propagation.
    # This is the main loop for Dijkstra's algorithm and the fast marching method.            
    def perform_loop(self, order, grid, start_pix, end_flag, refine_flag=False, main_flag=False):
        
        if end_flag == True:
            self.end_flag = True
        else:
            self.end_flag = False
        self.geo_distances[start_pix.pixel_index] = 0.0
        if main_flag == False:
            self.queue = [(0.0, start_pix.pixel_index)]
            no_alive = 0 # Number of pixels with a set distance
        else:
            no_alive = len([i for i in self.alive.values() if i == True])
        
        end_prev = -1
        end_neighbour = -1
        
        # Perform wavefront propagation
        while no_alive != grid.no_pixels:
            # Obtain index of pixel with smallest distance
            trial_distance, trial = heapq.heappop(self.queue)
            # If obtained trial index already alive, skip iteration
            if self.alive[trial] == True:
                continue
            
            # Set trial pixel to alive
            self.alive[trial] = True
            no_alive += 1
            
            # Check end point
            if self.end_flag == True:
                if trial in self.end_pixel_dict:
                    for k in self.end_pixel_dict[trial]:
                        self.no_endpoints_alive += 1
                    #self.end_distance = self.fmm_distance(grid, end_prev, self.end_pixel.pixel_index, end_neighbour, order, start_pix)
                    if self.no_endpoints_alive == len(self.end_pixel):
                        break
            # Check border (refined only)
           # if refine_flag == True and trial in grid.border_pixels:
            if refine_flag == True and grid.pixel[trial].is_border == True:
                break
            
            
            # Visit neighbours of trial pixel
            for i in range(len(grid.pixel[trial].neighbour)):
                visit = grid.pixel[trial].neighbour[i]
                if self.alive[visit] == True:
                    continue
                proposed_distance = self.fmm_distance(grid, trial, grid.pixel[visit], i, order, start_pix)
                if proposed_distance < self.geo_distances[visit]:
                    # Add visited pixel to the queue and update value in geo_distances
                    self.geo_distances[visit] = proposed_distance
                    heapq.heappush(self.queue, (proposed_distance, visit))
                    # If appropriate, record pixel used to derive this distance
                    if self.end_flag == True:
                        if visit in self.end_pixel_dict:
                            for k in self.end_pixel_dict[visit]:
                                self.end_distance[k] = self.fmm_distance(grid, trial, self.end_pixel[k], i, order, start_pix)
                                #end_prev = trial
                                #end_neighbour = i

        
        
    # Determine distance from a trial to visited pixel
    def fmm_distance(self, grid, trial, pix_visit, neighbour_no, order, start_pix):
        if order > 0: # Fast marching method
            if pix_visit.pixel_index == 0 or pix_visit.pixel_index == grid.no_pixels-1: # If visited pixel is a pole
                proposed_distance = -1.0
            else:
                #print(trial, ' , ', pix_visit)
                proposed_distance = self.fast_marching_method(pix_visit, order, grid, start_pix)
            # If the fast marching method is not applicable, get a distance of -1
            # and Dijkstra's algorithm (point-to-point distance) is used instead.
            if proposed_distance == -1.0:
                proposed_distance = self.point_to_point_distance(trial, pix_visit, neighbour_no, grid, start_pix)
        else:         # Dijkstra's algorithm
            proposed_distance = self.point_to_point_distance(trial, pix_visit, neighbour_no, grid, start_pix)
        return proposed_distance
    
    # Calculate distance between representative points of trial and neighbour
    def point_to_point_distance(self, trial, pix_visit, neighbour_no, grid, start_pix):
        if trial == start_pix.pixel_index:
            proposed_distance = start_pix.neighbour_distance[neighbour_no]
        else:
            if pix_visit.is_endpoint == True:
                #proposed_distance = self.get_neighbour_distance(trial, neighbour_no, grid, start_pix)
                proposed_distance = pix_visit.neighbour_distance[self.inverted_neighbour[neighbour_no]]
            else:
                proposed_distance = grid.get_distance(grid.pixel[trial], neighbour_no)
        '''    
            if self.end_flag == True:
                if grid.pixel[trial].neighbour[neighbour_no] == self.end_pixel.pixel_index:
                    proposed_distance = self.get_neighbour_distance(trial, neighbour_no, grid, start_pix)
                else:
                    proposed_distance = grid.get_distance(grid.pixel[trial], neighbour_no)
            else:
                proposed_distance = grid.get_distance(grid.pixel[trial], neighbour_no)
        '''
            
                
        proposed_distance += self.geo_distances[trial]
        
        return proposed_distance
    
    # Determine neighbour-to-neighbour distance (including start/end pixel cases)
    def get_neighbour_distance(self, pix_visit, neighbour_no, grid, start_pix):
        neighbour_dist = 0.0
        if pix_visit.neighbour[neighbour_no] == start_pix.pixel_index:
            # Get neighbour number of visit pixel relative to start pixel
            if pix_visit.neighbour[neighbour_no] == 0: # Neighbour is north pole
                neighbour_dist = start_pix.neighbour_distance[pix_visit.pixel_index-1]
            elif pix_visit.neighbour[neighbour_no] == grid.no_pixels-1: # Neighbour is south pole
                neighbour_dist = start_pix.neighbour_distance[grid.no_pixels-2-pix_visit.pixel_index]
            else: # Neighbour not a pole
                neighbour_dist = start_pix.neighbour_distance[inverted_neighbour[neighbour_no]]
        #elif self.end_flag == True:
        #    if grid.pixel[visit].neighbour[neighbour_no] == self.end_pixel.pixel_index:
        #        neighbour_dist = self.end_pixel.neighbour_distance[neighbour_no]
        #    else:
        #        neighbour_dist = grid.get_distance(grid.pixel[visit], neighbour_no)
        else:
            neighbour_dist = grid.get_distance(pix_visit, neighbour_no)
        return neighbour_dist
    
    
    # Calculate distances using the fast marching method
    def fast_marching_method(self, pix_visit, order, grid, start_pix):
        # Check for existence of 'alive' neighbours
        if ( ( self.alive[pix_visit.neighbour[0]] == True or self.alive[pix_visit.neighbour[1]] == True  ) and
             ( self.alive[pix_visit.neighbour[2]] == True or self.alive[pix_visit.neighbour[3]] == True  ) ):
            # i.e. If up or down neighbour alive AND left or right neighbour alive, use fast marching method
            
            # Determine whether Up or Down neighbour has smallest distance from start - use for theta distance.
            theta_index_1 = pix_visit.neighbour[0]
            neighbour_theta = 0     # Up
            dist_theta = self.geo_distances[theta_index_1]
            if self.geo_distances[pix_visit.neighbour[1]] < dist_theta:
                theta_index_1 = pix_visit.neighbour[1]
                dist_theta = self.geo_distances[theta_index_1]
                neighbour_theta = 1 # Down
            
            # Determine whether Left or Right neighbour has smallest distance from start - use for phi distance. 
            phi_index_1 = pix_visit.neighbour[2]
            neighbour_phi = 2     # Left
            dist_phi = self.geo_distances[phi_index_1]
            if self.geo_distances[pix_visit.neighbour[3]] < dist_phi:
                phi_index_1 = pix_visit.neighbour[3]
                dist_phi = self.geo_distances[phi_index_1]
                neighbour_phi = 3 # Right
            
            # Get distances from visited pixel to these neighbours
            delta_theta = self.get_neighbour_distance(pix_visit, neighbour_theta, grid, start_pix)
            delta_phi = self.get_neighbour_distance(pix_visit, neighbour_phi, grid, start_pix)
            
            # Use the fast marching method to calculate distance from start to visited pixel
            if order == 2:
                # 2nd Order FMM
                #
                # Check for second order theta neighbours
                pole_flag = True
                theta_index_2 = -1
                dist_theta_2 = dist_theta + 1
                if theta_index_1 != 0 and theta_index_1 != grid.no_pixels-1: # i.e. theta_index_1 is not a pole
                    theta_index_2 = grid.pixel[theta_index_1].neighbour[neighbour_theta]
                    dist_theta_2 = self.geo_distances[theta_index_2]
                    pole_flag = False
                
                # theta terms of quadratic
                if pole_flag == False and self.alive[theta_index_2] == True and dist_theta_2 <= dist_theta:
                    delta_theta_2 = self.get_neighbour_distance(grid.pixel[theta_index_1], neighbour_theta, grid, start_pix)
                    terms_theta = self.quadratic_terms(dist_theta, dist_theta_2, delta_theta, delta_theta_2, 2)
                else:
                    terms_theta = self.quadratic_terms(dist_theta, 0.0, delta_theta, 0.0, 1)
                    
                # Check for second order phi neighbours
                phi_index_2 = grid.pixel[phi_index_1].neighbour[neighbour_phi]
                dist_phi_2 = self.geo_distances[phi_index_2]
                
                # Phi terms of quadratic
                if self.alive[phi_index_2] and dist_phi_2 <= dist_phi:
                    delta_phi_2 = self.get_neighbour_distance(grid.pixel[phi_index_1], neighbour_phi, grid, start_pix)
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
                #
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
    
    # Find list of main grid pixels to use in refinement
    def find_refinement_pixels(self):
        # pix_refine is a list containing all main grid pixels within the 
        # refinement range of the start pixel.
        pix_refine = []
        if self.start_pixel.pixel_index == 0:                            # North pole
            pix_refine.append(0)
            for i in range(self.refine_range):
                for j in range(self.grid_main.n_phi):
                    pix_refine.append(self.grid_main.get_pixel_index(i+1, j))
        elif self.start_pixel.pixel_index == self.grid_main.no_pixels-1: # South pole
            pix_refine.append(self.grid_main.no_pixels-1)
            for i in range(self.refine_range):
                for j in range(self.grid_main.n_phi):
                    pix_refine.append(self.grid_main.get_pixel_index(self.grid_main.no_theta - 2 - i, j))
        else:                                                            # Non-pole
             for i in range(self.grid_main.no_pixels):
                 th = self.grid_main.pixel[i].theta_index
                 ph = self.grid_main.pixel[i].phi_index
                 if ( (abs(th - self.start_pixel.theta_index) <= self.refine_range) and
                    (  abs(ph - self.start_pixel.phi_index) <= self.refine_range or
                       abs(ph - self.start_pixel.phi_index) >= self.grid_main.no_phi - self.refine_range or
                       self.grid_main.pixel[i].is_north == True or self.grid_main.pixel[i].is_south == True) ):
                     pix_refine.append(i)  
        return pix_refine
            
            
    # Initialise a refined grid for use in a source refined fast marching method.
    def initialise_refined_grid(self, grid):
        
        # Create refined grid
        #no_theta_rfnd = grid.no_theta * self.refine_theta - 2
        no_theta_rfnd = self.refine_theta*(grid.no_theta-1) + 1
        no_phi_rfnd = grid.no_phi * self.refine_phi
        no_theta_border = int(self.refine_range*self.refine_theta + (self.refine_theta-1)/2) - 1
        no_phi_border = int(self.refine_range*self.refine_phi + (self.refine_phi-1)/2) - 1
        self.grid_rfnd = GeoGrid(self.grid_main.shape, no_theta_rfnd, no_phi_rfnd, True)
        # Find (theta, phi) coordinates of the centre of the main grid start pixel
        # This is used to define the centre of the refined region
        th_i = grid.pixel[self.start_pixel.pixel_index].theta_index
        ph_i = grid.pixel[self.start_pixel.pixel_index].phi_index
        centre_theta = grid.theta_list[th_i]
        centre_phi = grid.phi_list[ph_i]
        # Find indices of the centre pixel, in the refied grid
        centre_theta = self.grid_rfnd.find_theta_index(centre_theta)
        centre_phi = self.grid_rfnd.find_phi_index(centre_phi)
        
        # Initialise refined grid
        self.pix_refine = self.find_refinement_pixels()
        self.main2rfnd = [-1] * len(self.pix_refine)
        self.grid_rfnd.pixel = {}
        self.grid_rfnd.refined_indices = []
        self.alive = {}
        self.geo_distances = {}
        self.no_border_pixels = 0
        
        if len(self.pix_refine) == self.grid_main.no_pixels:
            # Include all refined grid pixels
            for i in range(self.grid_rfnd.no_pixels):
                self.alive[i] = False
                self.geo_distances[i] = math.inf
                th = self.grid_rfnd.get_theta_index(i)
                ph = self.grid_rfnd.get_phi_index(i, th)
                carts = self.grid_rfnd.polars_to_cartesians( self.grid_rfnd.theta_list[th],
                                                             self.grid_rfnd.phi_list[ph] )
                self.grid_rfnd.pixel[i] = GeoPixel(i, th, ph, self.grid_rfnd.no_pixels, carts, is_endpoint=False)
                self.grid_rfnd.find_neighbour_indices(i, th, ph)
                
                # Find the relevant pixels in the refined grid
                for k in range(len(self.pix_refine)):
                    th_in = self.grid_main.pixel[self.pix_refine[k]].theta_index
                    ph_in = self.grid_main.pixel[self.pix_refine[k]].phi_index
                    self.main2rfnd[k] = self.grid_rfnd.find_pixel_index(self.grid_main.theta_list[th_in],
                                                                        self.grid_main.phi_list[ph_in])
        else:
            # Find the relevant pixels in the refined grid
            for k in range(len(self.pix_refine)):

                # Get (main grid) indices and values for theta and phi
                th_in = self.grid_main.pixel[self.pix_refine[k]].theta_index
                ph_in = self.grid_main.pixel[self.pix_refine[k]].phi_index
                #th = self.grid_main.theta_list[th_in]
                #ph = self.grid_main.phi_list[ph_in]
                # main2rfnd maps main grid pixels corresponding to pix_refine to the
                # equivalent refined grid pixel.
                self.main2rfnd[k] = self.grid_rfnd.find_pixel_index(self.grid_main.theta_list[th_in],
                                                                    self.grid_main.phi_list[ph_in])
                # Get theta and phi indices of central pixel
                # Find limits for loops over theta and phi (based on refine_theta and refine_phi)
                th_c = self.grid_rfnd.get_theta_index(self.main2rfnd[k])
                th_dn = th_c - int((self.refine_theta - 1)/2)
                th_up = th_c + int((self.refine_theta - 1)/2)
                if th_c == 0 or th_c == self.grid_rfnd.no_theta-1:
                    ph_c = 0
                    ph_dn = ph_c - self.refine_range*self.refine_phi + int((self.refine_phi - 1)/2)
                    ph_up = ph_c + self.refine_range*self.refine_phi + int((self.refine_phi - 1)/2)
                else:
                    ph_c = self.grid_rfnd.get_phi_index(self.main2rfnd[k], th_c)
                    ph_dn = ph_c - int((self.refine_phi - 1)/2)
                    ph_up = ph_c + int((self.refine_phi - 1)/2)
                    
                for i in range(th_dn, th_up+1):
                    # Ignore if theta out of range
                    if i >= 0 and i < self.grid_rfnd.no_theta:
                        
                        if i == 0:
                            # North pole
                            carts = np.zeros(3)
                            carts[2] = self.grid_rfnd.shape.c_axis
                            self.grid_rfnd.pixel[0] = GeoPixel(0, 0, 0, self.grid_rfnd.no_pixels, carts, is_endpoint=False)
                            self.grid_rfnd.find_neighbour_indices(0, 0, 0)
                            if centre_theta != 0 and centre_theta < no_theta_border:
                                self.grid_rfnd.pixel[0].is_border = True
                                self.no_border_pixels += 1
                            self.grid_rfnd.refined_indices.append(0)
                            self.alive[0] = False
                            self.geo_distances[0] = math.inf
                                
                        elif i == self.grid_rfnd.no_theta - 1:
                            # South pole
                            carts = np.zeros(3)
                            carts[2] = -self.grid_rfnd.shape.c_axis
                            self.grid_rfnd.pixel[self.grid_rfnd.no_pixels - 1] = GeoPixel(self.grid_rfnd.no_pixels - 1, self.grid_rfnd.no_theta - 1, 0, self.grid_rfnd.no_pixels, carts, is_endpoint=False)
                            self.grid_rfnd.find_neighbour_indices(self.grid_rfnd.no_pixels - 1, self.grid_rfnd.no_theta - 1, 0)
                            if centre_theta != self.grid_rfnd.no_theta - 1 and self.grid_rfnd.no_theta - 1 - centre_theta < no_theta_border:
                                self.grid_rfnd.pixel[self.grid_rfnd.no_pixels - 1].is_border = True
                                self.no_border_pixels += 1
                            self.grid_rfnd.refined_indices.append(self.grid_rfnd.no_pixels - 1)
                            self.alive[self.grid_rfnd.no_pixels - 1] = False
                            self.geo_distances[self.grid_rfnd.no_pixels - 1] = math.inf
                            
                        else:
                            for j in range(ph_dn, ph_up+1):
                                
                                # Account for periodicity of phi
                                if i == 0 or i == self.grid_rfnd.no_theta-1:
                                    j_in = 0 # Pixel is a pole
                                elif j < 0:
                                    j_in = j + self.grid_rfnd.no_phi
                                elif j >= self.grid_rfnd.no_phi:
                                    j_in = j - self.grid_rfnd.no_phi
                                else:
                                    j_in = j
                                # Find index of pixel
                                if i == th_c and j == ph_c:
                                    # Central pixel
                                    pix_in = self.main2rfnd[k]
                                else:
                                    pix_in = self.grid_rfnd.get_pixel_index(i, j_in)
                                # Create pixel
                                if self.is_flat == False:
                                    carts = self.grid_rfnd.polars_to_cartesians( self.grid_rfnd.theta_list[i],
                                                                                 self.grid_rfnd.phi_list[j_in] )
                                else:
                                    carts = [self.grid_rfnd.phi_list[j_in], self.grid_rfnd.theta_list[i]]
                                
                                self.grid_rfnd.pixel[pix_in] = GeoPixel(pix_in, i, j_in, self.grid_rfnd.no_pixels, carts, is_endpoint=False)
                                self.grid_rfnd.find_neighbour_indices(pix_in, i, j_in)
                                
                                # Check whether pixel is on refined region border
                                th_diff = abs(i - centre_theta)
                                ph_diff = j_in - centre_phi
                                # Handle phi periodicity
                                if ph_diff > self.grid_rfnd.no_phi/2:
                                    ph_diff -= self.grid_rfnd.no_phi
                                elif ph_diff < -self.grid_rfnd.no_phi/2:
                                    ph_diff += self.grid_rfnd.no_phi
                                # Check for border
                                if  th_diff == no_theta_border or abs(ph_diff) == no_phi_border:
                                   self.grid_rfnd.pixel[pix_in].is_border = True
                                   self.no_border_pixels += 1
                                   
                                self.grid_rfnd.refined_indices.append(pix_in)
                                self.alive[pix_in] = False
                                self.geo_distances[pix_in] = math.inf
                        
            # Get number of refined pixels and initialise arrays with this length
            self.no_refined_pixels = len(self.grid_rfnd.refined_indices)
        
        # Create refined grid starting pixel
        self.start_pixel_rfnd = self.initialise_pixel(self.theta, self.phi, self.grid_rfnd) 
        
    # Go from refined grid to the main grid, ready to proceed with FMM calculation
    def transfer_grid(self):
        # Map refined grid data to the main grid
        
        # Initialise temporary dictionaries for main grid data
        alive_main = {}
        geo_distances_main = {}
        for i in range(self.grid_main.no_pixels):
            alive_main[i] = False
            geo_distances_main[i] = math.inf
        self.queue = []
        
        for i in range(len(self.pix_refine)):
            pix_main = self.pix_refine[i] # Index of a relevant main grid pixel
            pix_rfnd = self.main2rfnd[i]  # Index of corresponding refined grid pixel
            # Transfer distance from start point
            geo_distances_main[pix_main] = self.geo_distances[pix_rfnd]
            
            if self.alive[pix_rfnd] == True:
                count = 0 # Number of neighbours that are alive on the refined grid
                count_d = 0 # Number of neighbours with non-infinite distance
                for pix_nei in self.grid_main.pixel[pix_main].neighbour:
                    if pix_nei in self.pix_refine:
                        if self.alive[self.main2rfnd[self.pix_refine.index(pix_nei)]] == True:
                            count += 1
                        if self.geo_distances[self.main2rfnd[self.pix_refine.index(pix_nei)]] < math.inf:
                            count_d += 1

                if count == 4:
                    # Pixel is set to alive if all 4 neighbours have a set value
                    #alive_main.append(pix_main)
                    alive_main[pix_main] = True
                elif count_d > 1:
                    # Otherwise, add to the queue
                    heapq.heappush(self.queue, (geo_distances_main[pix_main], pix_main))
            #elif self.geo_distances[pix_rfnd] < math.inf:
                #heapq.heappush(self.queue, (geo_distances_main[pix_main], pix_main))

        # Copy 'main grid' alive and geo_distances arrays to the GeoFMM version
        
        c = 0
        for i in range(len(alive_main.keys())):
            if alive_main[i] == True:
                c += 1
        
        self.alive = copy.copy(alive_main)
        self.geo_distances = copy.copy(geo_distances_main)
                    
         
    # Get distance to pixel
    def get_distance(self, pixel_index):
        if self.is_refined_grid == True:
            return self.geo_distances[self.grid_rfnd.refined_indices.index(pixel_index)]
        else:
            return self.geo_distances[pixel_index]
        
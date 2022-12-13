# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:29:15 2022

@author: Callum Marples
"""

import math
import numpy as np
import heapq

from math import sqrt

# This class contains information for a vertex in the grid used in the fast
# marching method.
class FmmResult:
    def __init__(self, no_vertices):
        self.distance = [math.inf] * no_vertices
        self.accepted = [False] * no_vertices
        self.update = [[-1]] * no_vertices
        self.no_obtuse = 0
        
def fast_marching(grid, start_vertex, start_carts, is_dijkstra=False, end_dict={}):
    # Initialise output object
    fmm = FmmResult(grid.no_vertices)
    
    # Start point
    d_start = ( (grid.vertex[start_vertex].carts[0] - start_carts[0])**2.0 +
              (  grid.vertex[start_vertex].carts[1] - start_carts[1])**2.0 +
              (  grid.vertex[start_vertex].carts[2] - start_carts[2])**2.0 )
    if d_start < 1e-12: # Start point coincides with a vertex
        fmm.distance[start_vertex] = 0.0
        # Begin queue
        queue = [(0.0, start_vertex)]
    else: # Start point is not a vertex
        d_start = math.sqrt(d_start)
        fmm.distance[start_vertex] = d_start
        queue = [(d_start, start_vertex)]
        for j in grid.vertex[start_vertex].neighbour.keys():
            dist = euclidean_distance(grid.vertex[j].carts, start_carts)
            fmm.distance[j] = dist
            heapq.heappush(queue, (dist, j))
        
    
    no_accepted = 0

    no_ends = len(end_dict)
    no_ends_accepted = 0
    
    # Perform wavefront propagation
    while no_accepted != grid.no_vertices:
        
        # Obtain index of vertex with smallest distance from the start
        trial_distance, trial = heapq.heappop(queue)
        # If obtained trial index already accepted, skip iteration
        if fmm.accepted[trial] == True:
            continue
        # Add trial vertex to accepted set
        fmm.accepted[trial] = True
        no_accepted += 1
        if no_ends > 0 and trial in end_dict.keys():
            end_dict[trial] = True
            no_ends_accepted += 1
            if no_ends_accepted == no_ends:
                break
        
        # Visit neighbours of trial vertex
        for visit in grid.vertex[trial].neighbour_set:
            if fmm.accepted[visit] == True:
                continue
            proposed_distance = fast_marching_update(visit, trial, grid.vertex, fmm, is_dijkstra)
            if proposed_distance < fmm.distance[visit]:
                # Add visited vertex to the queue and update value in fmm.distance
                fmm.distance[visit] = proposed_distance
                heapq.heappush(queue, (proposed_distance, visit))
    return fmm

    
def fast_marching_update(visit, trial, vertex, fmm, is_dijkstra):
    if is_dijkstra == False:
        # Does an accepted third vertex exist?       
        support = get_supporting_vertex(visit, trial, vertex, fmm)
        
        if support == -1:
            fmm.update[visit] = [trial]
            #return get_distance(vertex, visit, trial) + fmm.distance[trial]
            return vertex[visit].neighbour[trial].distance + fmm.distance[trial]
        else:
            fmm.update[visit] = [trial, support]
            #v = fmm_first_order_update_general(visit, trial, support, vertex, vertex[visit].carts, fmm)
            v = fmm_first_order_update(visit, trial, support, vertex, fmm)
            if v[0] == -1:
                v[0] = vertex[visit].neighbour[trial].distance + fmm.distance[trial]
            return v[0]
    else:
        # Dijkstra's algorithm
        fmm.update[visit] = [trial]
        #return get_distance(vertex, visit, trial) + fmm.distance[trial]
        return vertex[visit].neighbour[trial].distance + fmm.distance[trial]
    
    
    
def get_supporting_vertex(visit, trial, vertex, fmm):
    v1 = vertex[visit].neighbour[trial].face[0]
    v2 = vertex[visit].neighbour[trial].face[1]
    # Are these vertices accepted?
    if fmm.accepted[v1] == True:
        if fmm.accepted[v2] == True:
            # Both vertices valid, so select the one closest to the start
            if fmm.distance[v1] < fmm.distance[v2]:
                support = v1
            else:
                support = v2
        else:
            # v1 accepted but not v2, so select v1
            support = v1
    else:
        if fmm.accepted[v2] == True:
            # v2 accepted but not v1, so select v2
            support = v2
        else:
            # Neither vertex accepted, so use point-to-point distance
            support = -1
    return support




    


def fmm_idw3(end_carts, carts_1, carts_2, carts_3, d1, d2, d3):
    w1 = 1.0 / ( (end_carts[0] - carts_1[0])**2.0 + (end_carts[1] - carts_1[1])**2.0 + (end_carts[2] - carts_1[2])**2.0 )
    w2 = 1.0 / ( (end_carts[0] - carts_2[0])**2.0 + (end_carts[1] - carts_2[1])**2.0 + (end_carts[2] - carts_2[2])**2.0 )
    w3 = 1.0 / ( (end_carts[0] - carts_3[0])**2.0 + (end_carts[1] - carts_3[1])**2.0 + (end_carts[2] - carts_3[2])**2.0 )
    return (w1*d1 + w2*d2 + w3*d3) / (w1 + w2 + w3)


def fmm_idw(end_carts, carts, dist):
    w = []
    num = 0.0
    den = 0.0
    for i in range(len(carts)):
        w.append( 1.0 / ( (end_carts[0] - carts[i][0])**2.0 + (end_carts[1] - carts[i][1])**2.0 + (end_carts[2] - carts[i][2])**2.0 ) )
        den += w[i]
        num += w[i] * dist[i]
    return num / den

# Use FMM data to obtain a distance to an endpoint
def endpoint_distance(vertex, fmm, end_carts, end_vertex, shape):

    # Find neighbour information and interpolate
    dist = [fmm.distance[end_vertex]]
    carts = [vertex[end_vertex].carts]
    for j in vertex[end_vertex].neighbour.keys():
        carts.append(vertex[j].carts)
        dist.append(fmm.distance[j])
    return fmm_idw(end_carts, carts, dist)
    
    

# Calculate Euclidean distance between two points x and y
def euclidean_distance(x, y):
    return math.sqrt( (x[0] - y[0])**2.0 + (x[1] - y[1])**2.0 + (x[2] - y[2])**2.0 )







def fmm_first_order_update(visit, trial, support, vertex, fmm):
    
    cos_psi = vertex[visit].neighbour[trial].face_angle[support]
    
    if cos_psi < 0.0:
        fmm.no_obtuse += 1
    
    if cos_psi < 0.0:
        
        fmm.update[visit] = [trial, -1]
        return [-1, 0, 0, 0]
    
    else:
        
        a = vertex[visit].neighbour[trial].distance
        b = vertex[visit].neighbour[support].distance
        
        u = fmm.distance[trial] - fmm.distance[support]
        
        asq = a*a
        bsq = b*b
        b2 = 2.0*b
        a_cos_psi = a*cos_psi
        alpha = asq + bsq - b2*a_cos_psi
        beta = b2*u*(a_cos_psi - b)
        gamma = bsq*(u*u - asq + a_cos_psi*a_cos_psi)
        
        # Solve quadratic
        discr = beta*beta - 4.0*alpha*gamma
        if discr > 0.0:
            t = ( -beta + sqrt(discr) ) / (2.0*alpha)
    
            # Check upwinding condition
            x = b*(t-u)/t
            if u < t and ((cos_psi >= 0.0 and cos_psi < 1.0e-12) or (x > a_cos_psi and x < a/cos_psi)):
                fmm.update[visit] = [trial, support]
                return [t + fmm.distance[support], alpha, beta, gamma]
            else:
                fmm.update[visit] = [trial]
                return [-1, alpha, beta, gamma]
        else:
            fmm.update[visit] = [trial]
            return [-1, alpha, beta, gamma]
    
    
    
    
    
    
    
    
    
def fmm_first_order_update_general(visit, trial, support, vertex, visit_carts, fmm):
    
    # Vectors of triangle sides
    p1 = (vertex[trial].carts - visit_carts)
    p2 = (vertex[support].carts - visit_carts)
    p11 = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]
    p12 = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]
    p22 = p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]
    
    if p12 < 0:
        print('obtuse')
    
    # Find vectors a and b
    a = np.array([1.0, 1.0])
    b = np.array([-fmm.distance[trial], -fmm.distance[support]])
    
    # Set up the quadratic A^2 u + Bu + C = 0
    #alpha = a[0]*(a[0]*p22 - a[1]*p12) + a[1]*(a[1]*p11 - a[0]*p12)
    #beta = 2.0 * (a[0]*(b[0]*p22 - b[1]*p12) + a[1]*(b[1]*p11 - b[0]*p12))
    alpha = p22 + p11 - 2.0*p12
    beta = 2.0 * ( b[0]*p22 + b[1]*p11 - p12*(b[0] + b[1]) )
    gamma = b[0]*(b[0]*p22 - b[1]*p12) + b[1]*(b[1]*p11 - b[0]*p12) - (p11*p22 - p12*p12)
    
    # Solve quadratic
    discr = beta*beta - 4.0*alpha*gamma
    if discr > 0:
        u = ( -beta + math.sqrt(discr) ) / (2.0*alpha)
        
        # Check upwinding condition  
        v1 = p22*(-b[0] - u) - p12*(-b[1] - u)
        v2 = p11*(-b[1] - u) - p12*(-b[0] - u)
        #v = math.sqrt(p11*p22)
        #v1 = v*(-b[0] - u) - p12*(-b[1] - u)
        #v2 = v*(-b[1] - u) - p12*(-b[0] - u)
        if u > -b[0] and u > -b[1] and v1 < 0 and v2 < 0:
            fmm.update[visit] = [trial, support]
            return [u, alpha, beta, gamma]
        else:
            fmm.update[visit] = [trial]
            return [-1, alpha, beta, gamma]
    else:
        fmm.update[visit] = [trial]
        return [-1, alpha, beta, gamma]
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:29:15 2022

@author: Callum Marples
"""

import math
import heapq

from math import sqrt

# This class contains information for a vertex in the mesh used in the fast
# marching method.
class FmmResult:
    def __init__(self, no_vertices):
        self.distance = [math.inf] * no_vertices
        self.accepted = [False] * no_vertices
        self.update = [[-1]] * no_vertices
        self.no_obtuse = 0
        
def fast_marching(mesh, start_vertex, start_carts, is_dijkstra=False, end_dict={}):
    # Initialise output object
    fmm = FmmResult(mesh.no_vertices)
    
    # Start point
    d_start = ( (mesh.vertex[start_vertex].carts[0] - start_carts[0])**2.0 +
              (  mesh.vertex[start_vertex].carts[1] - start_carts[1])**2.0 +
              (  mesh.vertex[start_vertex].carts[2] - start_carts[2])**2.0 )
    if d_start < 1e-12: # Start point coincides with a vertex
        fmm.distance[start_vertex] = 0.0
        # Begin queue
        queue = [(0.0, start_vertex)]
    else: # Start point is not a vertex
        d_start = math.sqrt(d_start)
        fmm.distance[start_vertex] = d_start
        queue = [(d_start, start_vertex)]
        for j in mesh.vertex[start_vertex].neighbour.keys():
            dist = euclidean_distance(mesh.vertex[j].carts, start_carts)
            fmm.distance[j] = dist
            heapq.heappush(queue, (dist, j))
        
    
    no_accepted = 0

    no_ends = len(end_dict)
    no_ends_accepted = 0
    
    # Perform wavefront propagation
    while no_accepted != mesh.no_vertices:
        
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
        for visit in mesh.vertex[trial].neighbour_set:
            if fmm.accepted[visit] == True:
                continue
            proposed_distance = fast_marching_update(visit, trial, mesh.vertex, fmm, is_dijkstra)
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
    
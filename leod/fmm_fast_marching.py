# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:29:15 2022

@author: Callum Marples
"""

import math
import numpy as np
import heapq

from leod.fmm_vertex import get_distance

# This class contains information for a vertex in the grid used in the fast
# marching method.
class FmmResult:
    def __init__(self, no_vertices):
        self.distance = [math.inf] * no_vertices
        self.accepted = [False] * no_vertices
        
def fast_marching(vertex, start, order):
    # Initialise output object
    no_vertices = len(vertex)
    fmm = FmmResult(no_vertices)
    # Start point
    fmm.distance[start] = 0.0
    # Begin queue
    queue = [(0.0, start)]
    no_accepted = 0
    
    # Perform wavefront propagation
    while no_accepted != no_vertices:
        
        # Obtain index of vertex with smallest distance from the start
        trial_distance, trial = heapq.heappop(queue)
        # If obtained trial index already accepted, skip iteration
        if fmm.accepted[trial] == True:
            continue
        # Add trial vertex to accepted set
        fmm.accepted[trial] = True
        no_accepted += 1
        
        # Visit neighbours of trial vertex
        for i in range(len(vertex[trial].neighbour)):
            visit = vertex[trial].neighbour[i]
            if fmm.accepted[visit] == True:
                continue
            proposed_distance = fast_marching_update(visit, trial, vertex, fmm, order)
            if proposed_distance < fmm.distance[visit]:
                # Add visited vertex to the queue and update value in fmm.distance
                fmm.distance[visit] = proposed_distance
                heapq.heappush(queue, (proposed_distance, visit))
    return fmm

    
def fast_marching_update(visit, trial, vertex, fmm, order):
    if order > 0:
        # Does an accepted third vertex exist?
        support = get_supporting_vertex(visit, trial, vertex, fmm)
        if support == -1:
            return get_distance(vertex, visit, trial) + fmm.distance[trial]
        else:
            v = fmm_first_order_update(visit, trial, support, vertex, fmm)
            return v[0]
    else:
        # Dijkstra's algorithm
        return get_distance(vertex, visit, trial) + fmm.distance[trial]
    
    
    
def get_supporting_vertex(visit, trial, vertex, fmm):
    v3 = [-1, -1]
    j = 0
    for i in range(len(vertex[visit].face)):
        if trial in vertex[visit].face[i]:
            if vertex[visit].face[i][0] == trial:
                v3[j] = vertex[visit].face[i][1]
            else:
                v3[j] = vertex[visit].face[i][0]
            j += 1
    v1 = v3[0]
    v2 = v3[1]
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



def fmm_first_order_update(visit, trial, support, vertex, fmm):
    
    # Vectors of triangle sides
    p1 = (vertex[trial].carts - vertex[visit].carts)
    p2 = (vertex[support].carts - vertex[visit].carts)
    p11 = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]
    p12 = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]
    p22 = p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]
    
    
    # Find vectors a and b
    a = np.array([1.0, 1.0])
    b = np.array([-fmm.distance[trial], -fmm.distance[support]])
    
    # Set up the quadratic A^2 u + Bu + C = 0
    alpha = a[0]*(a[0]*p22 - a[1]*p12) + a[1]*(a[1]*p11 - a[0]*p12)
    beta = 2.0 * (a[0]*(b[0]*p22 - b[1]*p12) + a[1]*(b[1]*p11 - b[0]*p12))
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
            return [u, alpha, beta, gamma]
        else:
            return [get_distance(vertex, visit, trial) + fmm.distance[trial], alpha, beta, gamma]
    else:
        return [get_distance(vertex, visit, trial) + fmm.distance[trial], alpha, beta, gamma]
    

def fmm_first_order_update_original(visit, trial, support, vertex, fmm):
    
    a = get_distance(vertex, visit, trial)
    b = get_distance(vertex, visit, support)
    w1 = (vertex[trial].carts - vertex[visit].carts) / a
    w2 = (vertex[support].carts - vertex[visit].carts) / b
    cos_psi = np.dot(w1, w2)
    u = fmm.distance[trial] - fmm.distance[support]
    
    alpha = a*a + b*b - 2*a*b*cos_psi
    beta = 2*b*u*(a*cos_psi - b)
    gamma = b*b*(u*u - a*a*(1-cos_psi*cos_psi))
    
    # Solve quadratic
    discr = beta*beta - 4.0*alpha*gamma
    if discr > 0.0:
        t = ( -beta + math.sqrt(discr) ) / (2.0*alpha)

        # Check upwinding condition
        x = b*(t-u)/t
        if u < t and x > a*cos_psi and x < a/cos_psi:
            return [t + fmm.distance[support], alpha, beta, gamma]
        else:
            return [get_distance(vertex, visit, trial) + fmm.distance[trial], alpha, beta, gamma]
    else:
        return [get_distance(vertex, visit, trial) + fmm.distance[trial], alpha, beta, gamma]
    
    
    
    
    
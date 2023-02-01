"""! 
@brief Routines and classes to perform the fast marching method.
@file fast_marching.py
@author Callum Marples

- Created on 11/10/2022. 
- Last modified on 31/01/2023.
"""

import heapq
from math import sqrt, inf, fabs

# This class contains information for a vertex in the mesh used in the fast
# marching method.
class FmmResult:
    """! The fast marching result class.
    Contains information obtained from applying the fast marching method to a FmmMesh object.
    """
    def __init__(self, no_vertices):
        self.distance = [inf] * no_vertices
        self.accepted = [False] * no_vertices
        self.update = [[-1]] * no_vertices
        self.no_obtuse = 0
        
def fast_marching(mesh, start_vertex, start_carts, is_dijkstra=False, end_dict={}):
    """! @brief Calculate approximate distances over a surface mesh, using either the fast marching method or Dijkstra's algorithm.
    @param mesh : FmmMesh \n
        The mesh on which the wavefront algorithm is applied.
    @param start_vertex : int \n
        Scalar index of the starting vertex.
    @param start_carts : 3-element NumPy array \n
        Cartesian coordinates of the start point. This need not coincide with the starting vertex.
    @param is_dijkstra : bool (optional)
        If True, Dijkstra's algorithm is used. Otherwise, the fast marching method is applied.
        Defaults to False.
    @param end_dict : dictionary of bool (optional)
        Dictionary containing any end points. The i-th entry corresponds to having endpoint(s) near vertex i.
        This allows the routine to terminate when the information for end points is found.
        Defaults to empty.
    @return FmmResult \n
        The results of the wavefront algorithm.
    """
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
        d_start = sqrt(d_start)
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



def euclidean_distance(x, y):
    """! @brief Compute the Euclidean distance between points \f$x\f$ and \f$y\f$.
    This is a subroutine of fast_marching.
    @param x : 3-element list of float \n
        Point \f$x\f$.
    @param y : 3-element list of float \n
        Point \f$y\f$.
    @return float \n
        Euclidean distance between \f$x\f$ and \f$y\f$.
    @see fast_marching
    """
    return sqrt( (x[0] - y[0])**2.0 + (x[1] - y[1])**2.0 + (x[2] - y[2])**2.0 )


    
def fast_marching_update(visit, trial, vertex, fmm, is_dijkstra):
    """! @brief Update the distance to a visited vertex. This is a subroutine of fast_marching.
    @param visit : int \n
        Index of the visited vertex (the vertex whose distance is updated).
    @param trial : int \n
        Index of the trial vertex (the vertex that has just been accepted).
    @param vertex : list of FmmVertex \n
        The list of vertices.
    @param fmm : FmmResult \n
        The object containing the results of the wavefront algorithm (this is updated within the routine).
    @param is_dijkstra : bool
        If True, Dijkstra's algorithm is used. Otherwise, the fast marching method is applied.
    @return float \n
        The updated distance from the start point to the visited vertex.
    @see fast_marching
    """
    if is_dijkstra == False:
        # Does an accepted third vertex exist?       
        support = get_supporting_vertex(visit, trial, vertex, fmm)
        
        if support == -1:
            fmm.update[visit] = [trial]
            return vertex[visit].neighbour[trial].distance + fmm.distance[trial]
        else:
            fmm.update[visit] = [trial, support]
            v = fmm_first_order_update(visit, trial, support, vertex, fmm)
            if v == -1:
                v = vertex[visit].neighbour[trial].distance + fmm.distance[trial]
            return v
    else:
        # Dijkstra's algorithm
        fmm.update[visit] = [trial]
        return vertex[visit].neighbour[trial].distance + fmm.distance[trial]
    
    
    
def get_supporting_vertex(visit, trial, vertex, fmm):
    """! @brief Given the visited and trial vertices, find the additional supporting
    vertex that completes the triangle for a fast marching update. This is a subroutine of fast_marching_update.
    @param visit : int \n
        Index of the visited vertex (the vertex whose distance is updated).
    @param trial : int \n
        Index of the trial vertex (the vertex that has just been accepted).
    @param vertex : list of FmmVertex \n
        The list of vertices.
    @param fmm : FmmResult \n
        The object containing the results of the wavefront algorithm (this is updated within the routine).
    @return int \n
        Index of the vertex used to complete the triangle for a fast marching update.
    @see fast_marching_update
    """
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



# Use FMM data to obtain a distance to an endpoint
def endpoint_distance(shape, vertex, fmm, end_carts, end_vertex):
    """! @brief Use fast marching data to obtain a distance to an end point.
    @param shape : EllipsoidShape \n
        The ellipsoid.
    @param vertex : list of FmmVertex \n
        The list of vertices.
    @param fmm : FmmResult \n
        The object containing the results of the wavefront algorithm.
    @param end_carts : 3-element NumPy array of float \n
        Cartesian coordinates of the end point.
    @param end_vertex : int \n
        Index of the vertex closest to the end point.
    @see distance_pair
    @see distance_multiple
    @see fmm_idw
    """
    # Find neighbour information and interpolate
    dist = [fmm.distance[end_vertex]]
    carts = [vertex[end_vertex].carts]
    for j in vertex[end_vertex].neighbour.keys():
        carts.append(vertex[j].carts)
        dist.append(fmm.distance[j])
    return fmm_idw(end_carts, carts, dist)



def fmm_idw(end_carts, carts, dist):
    """! @brief Interpolate the distance to an end point from fast marching data,
    using inverse square distance weighted interpolation.
    @param end_carts : 3-element NumPy array of float \n
        Cartesian coordinates of the end point.
    @param carts : list of 3-element NumPy array of float \n
        List of the Cartesian coordinates of vertices surrounding the end point.
    @param dist : list of float \n
        List of distances to the surrounding vertices obtained via the fast marching method.
    @return float
        Interpolated distance to the end point.
    @see endpoint_distance
    """
    w = []
    num = 0.0
    den = 0.0
    for i in range(len(carts)):
        w.append( 1.0 / ( (end_carts[0] - carts[i][0])**2.0 + (end_carts[1] - carts[i][1])**2.0 + (end_carts[2] - carts[i][2])**2.0 ) )
        den += w[i]
        num += w[i] * dist[i]
    return num / den



def fmm_first_order_update(visit, trial, support, vertex, fmm):
    """! @brief Given the visited and trial vertices, find the additional supporting
    vertex that completes the triangle for a fast marching update. This is a subroutine of fast_marching_update.
    @param visit : int \n
        Index of the visited vertex (the vertex whose distance is updated).
    @param trial : int \n
        Index of the trial vertex (the vertex that has just been accepted).
    @param support : int \n
        Index of the supporting vertex (the vertex that completes the triangle).
    @param vertex : list of FmmVertex \n
        The list of vertices.
    @param fmm : FmmResult \n
        The object containing the results of the wavefront algorithm (this is updated within the routine).
    @return float or int \n
        If float : The updated distance to the visited vertex. \n
        If int : A value of -1 is returned of the fast marching update violates monotonicity. In this case, Dijkstra's algorithm is used instead.
    @see fast_marching_update
    """
    cos_psi = vertex[visit].neighbour[trial].face_angle[support]
    
    if cos_psi < 0.0:
        fmm.no_obtuse += 1
    
    if cos_psi < 0.0:
        
        fmm.update[visit] = [trial, -1]
        return -1
    
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
            if fabs(t) > 1.0e-14:
                # Check upwinding condition
                x = b*(t-u)/t
                if u < t and ((cos_psi >= 0.0 and cos_psi < 1.0e-12) or (x > a_cos_psi and x < a/cos_psi)):
                    fmm.update[visit] = [trial, support]
                    return t + fmm.distance[support]
                else:
                    fmm.update[visit] = [trial]
                    return -1
            else:
                fmm.update[visit] = [trial]
                return -1
        else:
            fmm.update[visit] = [trial]
            return -1
    
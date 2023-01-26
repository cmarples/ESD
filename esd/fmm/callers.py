# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:38:56 2022

@author: Callum Marples

Prepare and call the fast marching routine for a given precaluculated graph and
input start/end points.
"""

from numpy import array
from numpy import zeros
from math import pi
from .mesh_ico import find_closest_vertex
from .fast_marching import fast_marching
from .fast_marching import endpoint_distance


def distance_pair(shape, mesh, start_point, end_point, is_dijkstra=False, is_radians=False):
    """! @brief Calculates the distance between given start and end points on an
                ellipsoid, with a given shape and mesh.
    @param shape : EllipsoidShape \n
        The ellipsoid.
    @param mesh : FmmMesh \n
        The mesh of vertices defined on the ellipsoid surface (can be polar or icosahedral).
    @param start_point : list of floats (2 elements) \n
        The start point, in \f$(\theta, \phi)\f$ coordinates.
    @param end_point : list of floats (2 elements) \n
        The end point, in \f$(\theta, \phi)\f$ coordinates.
    @param is_dijkstra : bool (optional) \n
        If True, Dijkstra's algorithm is used. Otherwise, the fast marching method
        is used. Defaults to False.
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return [float, FmmResult] \n
        The first element is the obtained distance from start_point to end_point.
        The second element contains more detailed information on the performed calculation.
    """
    start_point_temp = [0.0, 0.0]
    end_point_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
    for i in range(2):
        start_point_temp[i] = start_point[i] * conv
        end_point_temp[i] = end_point[i] * conv
    
    # Get Cartesian coordinatess of start and end points.
    start_carts = shape.polar2cart(start_point_temp[0], start_point_temp[1])   
    end_carts = shape.polar2cart(end_point_temp[0], end_point_temp[1])
     
    # Find start and end vertices
    if mesh.type == "pol":
        start_vertex = mesh.grid.find_vertex_index(start_point_temp[0], start_point_temp[1])
        end_vertex = mesh.grid.find_vertex_index(end_point_temp[0], end_point_temp[1])
    else:
        # Use triangulation
        start_vertex = find_closest_vertex(mesh, array(start_carts))
        end_vertex = find_closest_vertex(mesh, array(end_carts))
    
    
    end_vertex_carts = mesh.vertex[end_vertex].carts
    
    # If closest vertex and endpoint are the same, no interpolation needed.
    end_diff = end_vertex_carts - end_carts
    end_diff = end_diff[0]*end_diff[0] + end_diff[1]*end_diff[1] + end_diff[2]*end_diff[2]
    if end_diff > 1e-9:
        is_end_interpolate = True
    else:
        is_end_interpolate = False
    # Find neighbours of the endpoint, to use in the stopping criterion
    end_dict = {}
    end_dict[end_vertex] = False
    if is_end_interpolate == True:
        for j in mesh.vertex[end_vertex].neighbour.keys():
            end_dict[j] = False       
        
    # Call fast marching method
    fmm = fast_marching(mesh, start_vertex, start_carts, is_dijkstra, end_dict)
    
    if is_end_interpolate == True:
        d = endpoint_distance(shape, mesh.vertex, fmm, end_carts, end_vertex)
    else:
        d = fmm.distance[end_vertex]
    
    # Determine number of obtuse angles relevant to the endpoint distance
    temp_set = set()
    fmm.updater_set = set() # Set of vertices whose values influenced the endpoint
    for j in end_dict.keys():
        fmm.updater_set.add(j)
        temp_set.add(j)
    while len(temp_set) > 0:
        temp_set_new = set()
        for i in temp_set:
            for j in fmm.update[i]:
                if j != -1:
                    fmm.updater_set.add(j)
                    temp_set_new.add(j)
        temp_set = temp_set_new
    
    # Find obtuse vertices based on fmm.update
    fmm.no_obtuse_path = 0
    for i in fmm.updater_set:
        if len(fmm.update[i]) == 2 and fmm.update[i][1] == -1:
            fmm.no_obtuse_path += 1
    
    
    return d, fmm



def distance_multiple(shape, mesh, start_point, end_point=[], is_dijkstra=False, is_radians=False):
    """! @brief Calculates the distance between given start and end points on an
                ellipsoid, with a given shape and mesh.
    @param shape : EllipsoidShape \n
        The ellipsoid.
    @param mesh : FmmMesh \n
        The mesh of vertices defined on the ellipsoid surface (can be polar or icosahedral).
    @param start_point : list of floats (2 elements) \n
        The start point, in \f$(\theta, \phi)\f$ coordinates.
    @param end_point : \f$n\f$ element list of list of floats \n
        Set of \f$n\f$ end points, in \f$(\theta, \phi)\f$ coordinates.
    @param is_dijkstra : bool (optional) \n
        If True, Dijkstra's algorithm is used. Otherwise, the fast marching method
        is used. Defaults to False.
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return [list of floats, FmmResult] \n
        The first element is a list of obtained distances from the start_point to all points in end_point.
        The second element contains more detailed information on the performed calculation.
    """
    
    n = len(end_point)
        
    ### Convert to radians if input points given in degrees (assumed by default).
    start_point_temp = [0.0, 0.0]
    end_point_temp = zeros([n, 2])
    
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
        
    for j in range(2):
        start_point_temp[j] = start_point[j] * conv
    for i in range(n):
        for j in range(2):
            end_point_temp[i][j] = end_point[i][j] * conv
                
    ### Get Cartesian coordinates of start and end points.
    start_carts = shape.polar2cart(start_point_temp[0], start_point_temp[1]) 
    end_carts = [[0.0, 0.0, 0.0]] * n
    for i in range(n):
        end_carts[i] = shape.polar2cart(end_point_temp[i][0], end_point_temp[i][1])
    
    ### Find start and end vertices
    end_vertex = [-1] * n
    if mesh.type == "pol": # Polar mesh
        start_vertex = mesh.grid.find_vertex_index(start_point_temp[0], start_point_temp[1])
        for i in range(n):
            end_vertex[i] = mesh.grid.find_vertex_index(end_point_temp[i][0], end_point_temp[i][1])
    else: # Icosahedral triangulation
        start_vertex = find_closest_vertex(mesh, array(start_carts))
        for i in range(n):
            end_vertex[i] = find_closest_vertex(mesh, array(end_carts[i]))
   
    ### Prepare end points
    end_vertex_carts = [0] * n
    is_end_interpolate = [0] * n
    for i in range(n):
        end_vertex_carts[i] = mesh.vertex[end_vertex[i]].carts
        # If closest vertex and endpoint are the same, no interpolation needed.
        end_diff = end_vertex_carts[i] - end_carts[i]
        end_diff = end_diff[0]*end_diff[0] + end_diff[1]*end_diff[1] + end_diff[2]*end_diff[2]
        if end_diff > 1e-9:
            is_end_interpolate[i] = True
        else:
            is_end_interpolate[i] = False    
    
    # Find neighbours of the endpoints, to use in the stopping criterion
    end_dict = {}
    for i in range(n):
        end_dict[end_vertex[i]] = False
        if is_end_interpolate[i] == True:
            for j in mesh.vertex[end_vertex[i]].neighbour.keys():
                end_dict[j] = False  
            
    ### Call fast marching method
    fmm = fast_marching(mesh, start_vertex, start_carts, is_dijkstra, end_dict)
    
    ### Interpolate to find distance to endpoints
    d = [0.0] * n
    for i in range(n):
        if end_point[i][0] == start_point[0] and end_point[i][1] == start_point[1]:
            d[i] = 0.0
        else:
            if is_end_interpolate[i] == True:
                d[i] = endpoint_distance(shape, mesh.vertex, fmm, end_carts[i], end_vertex[i])
            else:
                d[i] = fmm.distance[end_vertex[i]]
            
    return d, fmm



def distance_all(shape, mesh, start_point, is_dijkstra=False, is_radians=False):
    """! @brief Calculates the distance between given a start point and all
    vertices on an ellipsoid, with a given shape and mesh.
    @param shape : EllipsoidShape \n
        The ellipsoid.
    @param mesh : FmmMesh \n
        The mesh of vertices defined on the ellipsoid surface (can be polar or icosahedral).
    @param start_point : list of floats (2 elements) \n
        The start point, in \f$(\theta, \phi)\f$ coordinates.
    @param is_dijkstra : bool (optional) \n
        If True, Dijkstra's algorithm is used. Otherwise, the fast marching method
        is used. Defaults to False.
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return [list of floats, FmmResult] \n
        The first element is a list of obtained distances from the start_point to all points in end_point.
        The second element contains more detailed information on the performed calculation.
    """
        
    ### Convert to radians if input points given in degrees (assumed by default).
    start_point_temp = [0.0, 0.0]
    
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
        
    for j in range(2):
        start_point_temp[j] = start_point[j] * conv
                
    ### Get Cartesian coordinates of start point.
    start_carts = shape.polar2cart(start_point_temp[0], start_point_temp[1]) 
    
    ### Find start vertex
    if mesh.type == "pol": # Polar mesh
        start_vertex = mesh.grid.find_vertex_index(start_point_temp[0], start_point_temp[1])
    else: # Icosahedral triangulation
        start_vertex = find_closest_vertex(mesh, array(start_carts)) 
            
    ### Call fast marching method
    fmm = fast_marching(mesh, start_vertex, start_carts, is_dijkstra, end_dict={})
    
    return fmm



def distance_end(shape, mesh, fmm, end_point, is_radians=False):
    """! @brief Interpolates the distance to a given end point on an ellipsoid 
    for which the FmmResult object has already been generated via distance_all.
    @param shape : EllipsoidShape \n
        The ellipsoid.
    @param mesh : FmmMesh \n
        The mesh of vertices defined on the ellipsoid surface (can be polar or icosahedral).
    @param fmm : FmmResult \n
        An output of the distance_all routine.
    @param end_point : list of floats (2 elements) \n
        The end point, in \f$(\theta, \phi)\f$ coordinates.
    @param is_radians : bool (optional) \n
        Specifies whether the elements start_point and end_point are given in
        radians (True) or degrees (False). If False, a coordiante conversion to
        radians is performed. Defaults to False.
    @return float \n
        Interpolated distance to the specified end point.
    @see distance_all
    """
    end_point_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = pi / 180.0
    else:
        conv = 1.0
    for i in range(2):
        end_point_temp[i] = end_point[i] * conv
    
    # Get Cartesian coordinatess of end point.  
    end_carts = shape.polar2cart(end_point_temp[0], end_point_temp[1])
     
    # Find end vertex
    if mesh.type == "pol":
        end_vertex = mesh.grid.find_vertex_index(end_point_temp[0], end_point_temp[1])
    else:
        # Use triangulation
        end_vertex = find_closest_vertex(mesh, array(end_carts))
    
    end_vertex_carts = mesh.vertex[end_vertex].carts
    
    # If closest vertex and endpoint are the same, no interpolation needed.
    end_diff = end_vertex_carts - end_carts
    end_diff = end_diff[0]*end_diff[0] + end_diff[1]*end_diff[1] + end_diff[2]*end_diff[2]
    if end_diff > 1e-9:
        is_end_interpolate = True
    else:
        is_end_interpolate = False
    # Find neighbours of the endpoint, to use in the stopping criterion
    end_dict = {}
    end_dict[end_vertex] = False
    if is_end_interpolate == True:
        for j in mesh.vertex[end_vertex].neighbour.keys():
            end_dict[j] = False       
    
    if is_end_interpolate == True:
        d = endpoint_distance(shape, mesh.vertex, fmm, end_carts, end_vertex)
    else:
        d = fmm.distance[end_vertex]
    
    return d
    
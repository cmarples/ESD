"""! 
@brief Icosahedral mesh routines.
@file mesh_ico.py
@author Callum Marples
- Created by Callum Marples on 27/10/2022.
- Last modified on 02/02/2022.
"""

import numpy as np
import random
from math import sqrt, fabs
from .classes import FmmVertex
from .classes import FmmNeighbour
from .classes import FmmMesh
from ..shape import EllipsoidShape

# Sphere triangulation using a geodesic polyhedron
def gen_ico_mesh(n=50, shape=EllipsoidShape(), is_split=False):
    """! @brief Generate a mesh approximating the surface of an ellipsoid, using the vertices of a regular icosahedron.
    @param n : int (optional) \n
        The number of points to take along each icosahedral edge. This parameter
        measures the resoltion of the mesh. Defaults to 50.
    @param shape : EllipsoidShape (optional) \n
        The ellipsoid. Defaults to an EllipsoidShape representing the unit sphere.
    @param is_split : bool (optional) \n
        If True, split each obtuse angle by adding an additional edge.
        Defaults to False.
    @return FmmMesh \n
        The mesh representing the surface of the ellipsoid specified by the input shape.
    """
    # Create FmmGrid object (includes empty vertex list).
    mesh = FmmMesh()
    
    # Icosahedron coordinates
    t = 0.5 * (1.0 + sqrt(5.0)) # Golden ratio
    
    ico_vertex = [ [0, 1, t], [0, -1, t], [0, 1, -t], [0, -1, -t],
                   [1, t, 0], [-1, t, 0], [1, -t, 0], [-1, -t, 0],
                   [t, 0, 1], [-t, 0, 1], [t, 0, -1], [-t, 0, -1] ]
    ico_vertex = np.array(ico_vertex)
    
    # Determine neighbouring vertices (edge length = 2)
    neigh = []
    for i in range(12):
        neigh.append([-1] * 5)
        neighbour_no = 0
        for j in range(12):
            if i != j:
                d2 = (ico_vertex[i][0]-ico_vertex[j][0])**2 + (ico_vertex[i][1]-ico_vertex[j][1])**2 + (ico_vertex[i][2]-ico_vertex[j][2])**2
                if fabs(d2 - 4.0) < 1e-9:
                    neigh[i][neighbour_no] = j
                    neighbour_no += 1
    
    # Determine faces by considering each trio in turn.
    ico_face = []
    for i in range(10):
        for j in range(i, 11):
            for k in range(j, 12):
                if j in neigh[i] and k in neigh[i] and i in neigh[j] and k in neigh[j] and i in neigh[k] and j in neigh[k]:
                    ico_face.append([i, j, k])
    
    # Determine the five faces surrounding each vertex.
    ico_vertex_face = []
    for j in range(12):
        ico_vertex_face.append([])
        for i in range(20):
            if j in ico_face[i]:
                ico_vertex_face[j].append(i)
                
    # Determine the two faces either side of each edge.
    ico_edge_face = {}
    for i in range(12):
        for j in neigh[i]:
            if j > i:
                ico_edge_face[(i, j)] = []
                for k in range(20):
                    if i in ico_face[k]:
                        if j in ico_face[k]:
                            ico_edge_face[(i, j)].append(k)
    
    # vertex_no : icosahedral vertex index
    def process_corner_vertex(vertex_no, index, ico_vertex_face):
        if ico_index[vertex_no] == -1:
            mesh.vertex.append(FmmVertex(index, ico_vertex[vertex_no]))
            mesh.vertex[-1].face = []
            mesh.vertex[-1].ico_face = []
            if j == 0:
                vertex_map[(j, 0)] = index
            else:
                vertex_map[(j, k)] = index
            ico_index[vertex_no] = index
            mesh.vertex[index].ico = ico_vertex_face[vertex_no] # Record icosahedral face that this vertex belongs to.
            index += 1
            
        else:
            if j == 0:
                vertex_map[(j, 0)] = ico_index[vertex_no]
            else:
                vertex_map[(j, k)] = ico_index[vertex_no]
        
        return index
    
    # vx, vy : icosahedral vertex indices
    def process_edge_vertex(vx, vy, jk, index, ico_edge_face):
        if vx < vy:
            x_lt_y = True
            edge = ico_edge_face[(vx, vy)]
            vertex_index = ico_edge_index[(vx, vy)][jk-1]
        else:
            x_lt_y = False
            edge = ico_edge_face[(vy, vx)]
            vertex_index = ico_edge_index[(vx, vy)][n-1-jk]                   
        if vertex_index == -1:
            mesh.vertex.append(FmmVertex(index, (v01*(j-k) + v02*k) / j))
            mesh.vertex[-1].face = []
            mesh.vertex[-1].ico_face = []
            vertex_map[(j, k)] = index
            if x_lt_y == True:
                ico_edge_index[(vx, vy)][jk-1] = index
            else:
                ico_edge_index[(vx, vy)][n-1-jk] = index
            
            mesh.vertex[index].ico = [edge[0], edge[1]] # Record icosahedral face that this vertex belongs to.
            index += 1 
            
        else:
            vertex_map[(j, k)] = vertex_index
        return index
                                
    ico_index = [-1] * 12 # Map from icosahedral index to triangulation index
    ico_edge_index = {}   # Map from an edge of a icosahedral face to an index
    for i in range(12):
        for j in neigh[i]:
            if j > i:
                ico_edge_index[(i, j)] = [-1] * (n-1)
    index = 0
    
    # Determine vertices, faces and neighbours of the triangulation
    for i in range(20):
        
        # Get icosahedron vertices
        v0 = ico_vertex[ico_face[i][0]]
        v1 = ico_vertex[ico_face[i][1]]
        v2 = ico_vertex[ico_face[i][2]]

        # Triangle vertices
        vertex_map = {}
        
        for j in range(n+1):
            
            v01 = (v0*(n-j) + v1*(j)) / n
            v02 = (v0*(n-j) + v2*(j)) / n
            
            if j == 0: # Vertex v0 may have already been defined.
                index = process_corner_vertex(ico_face[i][0], index, ico_vertex_face)
            else:
                for k in range(j+1):
                    if j == n:
                        if k == 0:
                            # Vertex v1 may have already been defined.
                            index = process_corner_vertex(ico_face[i][1], index, ico_vertex_face)
                        elif k == j:
                            # Vertex v2 may have already been defined.
                            index = process_corner_vertex(ico_face[i][2], index, ico_vertex_face)
                        else:
                            # Vertex along v1-v2 line may have already been defined
                            index = process_edge_vertex(ico_face[i][1], ico_face[i][2], k, index, ico_edge_face)
                    else:
                        if k == 0:
                            # Vertex along v0-v1 line may have already been defined
                            index = process_edge_vertex(ico_face[i][0], ico_face[i][1], j, index, ico_edge_face)    
                        elif k == j:
                            # Vertex along v0-v2 line may have already been defined
                            index = process_edge_vertex(ico_face[i][0], ico_face[i][2], j, index, ico_edge_face)
                        else:
                            # Vertex is within the icosahedral face
                            mesh.vertex.append(FmmVertex(index, (v01*(j-k) + v02*k) / j))
                            mesh.vertex[-1].face = []
                            mesh.vertex[-1].ico_face = []
                            vertex_map[(j, k)] = index
                            mesh.vertex[index].ico = [i] # Record icosahedral face that this vertex belongs to.
                            index += 1
                    
        # Find neighbours and faces for each vertex
        for j in range(n+1):
            if j == 0: # Vertex v0
                v = vertex_map[(j, 0)]
                mesh.vertex[v].neighbour[vertex_map[(1, 0)]] = FmmNeighbour()
                mesh.vertex[v].neighbour[vertex_map[(1, 1)]] = FmmNeighbour()
                mesh.vertex[v].face.append([vertex_map[(1, 1)], vertex_map[(1, 0)]])
            else:
                for k in range(j+1):
                    v = vertex_map[(j, k)]
                    if j == n:
                        if k == 0: # Vertex v1 
                            mesh.vertex[v].neighbour[vertex_map[(n-1, 0)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(n, 1)]] = FmmNeighbour()
                            mesh.vertex[v].face.append([vertex_map[(n-1, 0)], vertex_map[(n, 1)]])
                        elif k == j: # Vertex v2
                            mesh.vertex[v].neighbour[vertex_map[(n-1, n-1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(n, n-1)]] = FmmNeighbour()
                            mesh.vertex[v].face.append([vertex_map[(n-1, n-1)], vertex_map[(n, n-1)]])
                        else:
                            # Vertex along v1-v2 line
                            mesh.vertex[v].neighbour[vertex_map[(n-1, k)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(n, k+1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(n, k-1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(n-1, k-1)]] = FmmNeighbour()
                            mesh.vertex[v].face.append([vertex_map[(n-1, k)], vertex_map[(n, k+1)]])
                            mesh.vertex[v].face.append([vertex_map[(n, k-1)], vertex_map[(n-1, k-1)]])
                            mesh.vertex[v].face.append([vertex_map[(n-1, k-1)], vertex_map[(n-1, k)]])
                    else:
                        if k == 0:
                            # Vertex along v0-v1 line
                            mesh.vertex[v].neighbour[vertex_map[(j-1, 0)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j, 1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j+1, 1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j+1, 0)]] = FmmNeighbour()
                            mesh.vertex[v].face.append([vertex_map[(j-1, 0)], vertex_map[(j, 1)]])
                            mesh.vertex[v].face.append([vertex_map[(j, 1)], vertex_map[(j+1, 1)]])
                            mesh.vertex[v].face.append([vertex_map[(j+1, 1)], vertex_map[(j+1, 0)]])
                        elif k == j:
                            # Vertex along v0-v2 line
                            mesh.vertex[v].neighbour[vertex_map[(j+1, j+1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j+1, j)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j, j-1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j-1, j-1)]] = FmmNeighbour()
                            mesh.vertex[v].face.append([vertex_map[(j+1, j+1)], vertex_map[(j+1, j)]])
                            mesh.vertex[v].face.append([vertex_map[(j+1, j)], vertex_map[(j, j-1)]])
                            mesh.vertex[v].face.append([vertex_map[(j, j-1)], vertex_map[(j-1, j-1)]])
                        else:
                            # Vertex is within the icosahedral face
                            mesh.vertex[v].neighbour[vertex_map[(j-1, k)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j, k+1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j+1, k+1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j+1, k)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j, k-1)]] = FmmNeighbour()
                            mesh.vertex[v].neighbour[vertex_map[(j-1, k-1)]] = FmmNeighbour()
                            mesh.vertex[v].face.append([vertex_map[(j-1, k)], vertex_map[(j, k+1)]])
                            mesh.vertex[v].face.append([vertex_map[(j, k+1)], vertex_map[(j+1, k+1)]])
                            mesh.vertex[v].face.append([vertex_map[(j+1, k+1)], vertex_map[(j+1, k)]])
                            mesh.vertex[v].face.append([vertex_map[(j+1, k)], vertex_map[(j, k-1)]])
                            mesh.vertex[v].face.append([vertex_map[(j, k-1)], vertex_map[(j-1, k-1)]])
                            mesh.vertex[v].face.append([vertex_map[(j-1, k-1)], vertex_map[(j-1, k)]])
    
    # Number of vertices
    mesh.no_vertices = len(mesh.vertex)
    
    # Find possible faces for each neighbouring pair
    cos_alpha = 2.0 # Initialisation for face angles
    for i in range(mesh.no_vertices):
        for j in mesh.vertex[i].neighbour.keys():
            k_list = []
            for face_i in mesh.vertex[i].face:
                if j in face_i:
                    if face_i[0] == j:
                        k_list.append(face_i[1])                        
                    else:
                        k_list.append(face_i[0])
            
            mesh.vertex[i].neighbour[j].face = k_list
            mesh.vertex[i].neighbour[j].face_angle[k_list[0]] = cos_alpha
            mesh.vertex[i].neighbour[j].face_angle[k_list[1]] = cos_alpha
                    
    # Define outgoing neighbours
    for i in range(mesh.no_vertices):
        for j in mesh.vertex[i].neighbour.keys():
            mesh.vertex[i].neighbour_set.add(j)                        
                            
    # Project onto unit sphere surface
    for i in range(mesh.no_vertices):
        mesh.vertex[i].carts *= 1.0 / np.linalg.norm(mesh.vertex[i].carts)
    
    # Scale to ellipsoid
    for k in range(mesh.no_vertices):
        mesh.vertex[k].carts[0] *= shape.a_axis
        mesh.vertex[k].carts[1] *= shape.b_axis
        mesh.vertex[k].carts[2] *= shape.c_axis
        
    # Precalculate edge lengths and face angles
    mesh.precalculate()    
    
    if is_split == True:
        split_obtuse_angles(mesh)
    
    # Set grid type
    mesh.type = "ico"
    
    return mesh

# Find closest vertex to input point
def find_closest_vertex(mesh, p, v=0):
    """! @brief Find the closest vertex (in terms of Euclidean distance) to a given point.
    @param mesh : FmmMesh object \n
        The icosahedral mesh.
    @param p : 3-element NumPy array \n
        The point for which the closest vertex is to be found.
    @param v : int (optional)
        The starting vertex for the search. Defaults to 0.
    @return int \n
        The index of the vertex closest to point p.
    @see gen_ico_mesh
    @see distance_pair
    @see calculate_distances
    """
    d2 = {}
    is_min_found = False
    
    def perform_loop(vertex, p, v, d2, is_min_found):
        while is_min_found == False:
            if v not in d2.keys():
                vec = vertex[v].carts - p
                d2[v] = np.dot(vec, vec)
            u = v # Index of locally closest vertex
            d_min = d2[v]
            for j in vertex[v].neighbour.keys():
                if j not in d2.keys():
                    vec = vertex[j].carts - p
                    d2[j] = np.dot(vec, vec)
                if d2[j] < d_min:
                    d_min = d2[j]
                    u = j
            if u == v:
                is_min_found = True
            else:
                v = u
        return v, d_min
    
    max_dist = mesh.max_edge*mesh.max_edge # Largest neighbour distance squared.
    d_min = 2.0*max_dist
    while d_min > max_dist:
        v, d_min = perform_loop(mesh.vertex, p, v, d2, is_min_found)
        if d_min < max_dist:
            break
        else:
            # Select a random vertex for which d2 has not been calculated.
            while v in d2.keys():
                v = random.randint(0, len(mesh.vertex)-1)
    return v

def split_obtuse_angles(mesh):
    """! @brief Split each obtuse angle in a triangular mesh by adding additional edges.
    This routine updates the input vertex list and is a subroutine of gen_ico_mesh.
    @param mesh : FmmMesh object \n
        The icosahedral mesh.
    @see gen_ico_mesh
    """
    for i in range(len(mesh.vertex)):
        # Find faces to loop over
        face_set = set()
        for j in mesh.vertex[i].neighbour.keys():
            for k in mesh.vertex[i].neighbour[j].face:
                if k > j:
                    face_set.add((j, k))
        # Loop over each face, splitting if necessary
        for face in face_set:
            j = face[0]
            k = face[1]
            # Check the cosine of the angle to see if it is obtuse.
            if mesh.vertex[i].neighbour[j].face_angle[k] < 0.0:
                # Perform split on obtuse angle.
                # Find vertex s to perform split with.
                temp_set = mesh.vertex[j].neighbour.keys() & mesh.vertex[k].neighbour.keys()
                temp_set.remove(i)
                (s,) = temp_set
                # Add new neighbour
                mesh.vertex[i].neighbour[s] = FmmNeighbour()
                # Add outgoing neighbour
                mesh.vertex[s].neighbour_set.add(i)
                # Calculate distance
                mesh.vertex[i].neighbour[s].distance = sqrt( (mesh.vertex[i].carts[0]-mesh.vertex[s].carts[0])**2.0 +
                                                        (mesh.vertex[i].carts[1]-mesh.vertex[s].carts[1])**2.0 +
                                                        (mesh.vertex[i].carts[2]-mesh.vertex[s].carts[2])**2.0 )
                if mesh.vertex[i].neighbour[s].distance < mesh.min_edge:
                    mesh.min_edge = mesh.vertex[i].neighbour[s].distance
                if mesh.vertex[i].neighbour[s].distance > mesh.max_edge:
                    mesh.mxa_edge = mesh.vertex[i].neighbour[s].distance
                # Update face completion lists
                for v in range(2):
                    if mesh.vertex[i].neighbour[j].face[v] == k:
                        mesh.vertex[i].neighbour[j].face[v] = s
                    if mesh.vertex[i].neighbour[k].face[v] == j:
                        mesh.vertex[i].neighbour[k].face[v] = s
                mesh.vertex[i].neighbour[s].face = [j, k]
                # Remove old face angles        
                del mesh.vertex[i].neighbour[j].face_angle[k]
                del mesh.vertex[i].neighbour[k].face_angle[j]
                mesh.no_obtuse -= 1
                # Calculate new face angles
                wj = mesh.vertex[j].carts - mesh.vertex[i].carts
                wk = mesh.vertex[k].carts - mesh.vertex[i].carts
                ws = mesh.vertex[s].carts - mesh.vertex[i].carts
                dj = mesh.vertex[i].neighbour[j].distance
                dk = mesh.vertex[i].neighbour[k].distance
                ds = mesh.vertex[i].neighbour[s].distance
                mesh.vertex[i].neighbour[j].face_angle[s] = (wj[0]*ws[0] + wj[1]*ws[1] + wj[2]*ws[2]) / (dj*ds)
                mesh.vertex[i].neighbour[s].face_angle[j] = mesh.vertex[i].neighbour[j].face_angle[s]
                if mesh.vertex[i].neighbour[j].face_angle[s] < 0.0:
                    mesh.no_obtuse += 1
                mesh.vertex[i].neighbour[k].face_angle[s] = (wk[0]*ws[0] + wk[1]*ws[1] + wk[2]*ws[2]) / (dk*ds)
                mesh.vertex[i].neighbour[s].face_angle[k] = mesh.vertex[i].neighbour[k].face_angle[s]
                if mesh.vertex[i].neighbour[k].face_angle[s] < 0.0:
                    mesh.no_obtuse += 1
                
    # Determine the updated maximum and minimum face angles
    mesh.min_angle = -2.0
    mesh.max_angle = 2.0
    for i in range(mesh.no_vertices):        
            for j in mesh.vertex[i].neighbour:
                for k in mesh.vertex[i].neighbour[j].face_angle.keys():
                    if mesh.vertex[i].neighbour[j].face_angle[k] < mesh.max_angle:
                        mesh.max_angle = mesh.vertex[i].neighbour[j].face_angle[k]
                    if mesh.vertex[i].neighbour[j].face_angle[k] > mesh.min_angle:
                        mesh.min_angle = mesh.vertex[i].neighbour[j].face_angle[k]
                
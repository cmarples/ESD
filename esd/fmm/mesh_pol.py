# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:49:51 2022

@author: Callum Marples

Generate FmmGrid using an ellipsoid with scaled spherical polar coordinates
"""

import numpy as np

from .classes import FmmVertex
from .classes import FmmNeighbour
from .classes import FmmMesh
from ..shape import EllipsoidShape
from ..grid import Grid

def gen_pol_mesh(no_theta=91, no_phi=180, shape=EllipsoidShape(), is_connect_8=True):
    """! @brief Generate a mesh approximating the surface of an ellipsoid, using polar coordinates.
    @param no_theta : int (optional) \n
        The number of \f$\theta\f$ used in the mesh. Defaults to 91.
    @param no_phi : int (optional) \n
        The number of \f$\phi\f$ used in the mesh. Defaults to 180.
    @param shape : EllipsoidShape (optional) \n
        The ellipsoid. Defaults to an EllipsoidShape representing the unit sphere.
    @param is_connect_8 : bool (optional) \n
        If True, give each mesh vertex 8 neighbours (if applicable). If False, 4 neighbours are used.
        Defaults to True.
    @return FmmMesh \n
        The mesh representing the surface of the ellipsoid specified by the input shape.
    """
    # Create FmmMesh object (includes empty vertex list).
    mesh = FmmMesh()
    
    sphere = EllipsoidShape() # Unit sphere shape object.
    
    # Generate Grid object to store polar information.
    mesh.grid = Grid(no_theta, no_phi)
    mesh.no_vertices = (no_theta - 2)*no_phi + 2 # Used throughout FMM.
    
    # Calculate and store polar mesh information.
    mesh.no_theta = no_theta
    mesh.no_phi = no_phi
    
    # Construct the vertex array in the FmmGrid
    polar_index = []
    for i in range(mesh.no_vertices):
        
        # theta and phi indices
        th_index = mesh.grid.get_theta_index(i)
        ph_index = mesh.grid.get_phi_index(i, th_index)
        polar_index.append([th_index, ph_index])
        
        # Cartesian coordinates
        carts = np.array( sphere.polar2cart(mesh.grid.theta_list[th_index], mesh.grid.phi_list[ph_index]) )
        
        # Create new vertex
        mesh.vertex.append(FmmVertex(i, carts))
        
    # Find neighbouring vertices and update the graph
    for i in range(mesh.no_vertices):    
        find_neighbour_indices(mesh.vertex, i, polar_index[i][0], polar_index[i][1], mesh.grid, is_connect_8)
        
    # Define outgoing neighbours
    for i in range(mesh.no_vertices):
        for j in mesh.vertex[i].neighbour.keys():
            mesh.vertex[i].neighbour_set.add(j)   
    
    # Scale to ellipsoid
    for i in range(mesh.no_vertices):
        mesh.vertex[i].carts[0] *= shape.a_axis
        mesh.vertex[i].carts[1] *= shape.b_axis
        mesh.vertex[i].carts[2] *= shape.c_axis
        
    # Precalculate edge lengths and face angles
    mesh.precalculate()        
    
    # Set mesh type
    mesh.type = "pol"
    
    return mesh

   
# Find indices of neighbouring vertices
# Each pixel has four neighbours. For non-polar vertices, neighbours are
# always ordered as "Up, Down, Left, Right".
def find_neighbour_indices(vertex, i, th, ph, grid, is_connect_8=True):
    """! @brief Find and add neighbours for a given vertex. This is a subroutine of gen_pol_mesh.
    @param vertex : FmmVertex \n
        The vertex for which neighbours are to be found.
    @param i : int \n
        The scalar index of the vertex.
    @param th : int \n
        The \f$\theta\f$ index of the vertex.
    @param ph : int \n
        The \f$\phi\f$ index of the vertex.
    @param grid : Grid \n
        The \f$\theta\f$-\f$\phi\f$ grid used to generate the mesh.
    @param is_connect_8 : bool (optional) \n
        If True, give each mesh vertex 8 neighbours (if applicable). If False, 4 neighbours are used.
        Defaults to True.
    @see gen_pol_mesh
    """
    k0 = grid.no_vertices-1-grid.no_phi # First vertex at south pole adjacent band  
    cos_alpha = 2.0 # Initialisation for face angles
    
    # Pole neighbours
    if i == 0:
        
        # North pole
        vertex[0].neighbour[1] = FmmNeighbour()
        vertex[0].neighbour[1].face = [grid.no_phi, 2]
        vertex[0].neighbour[1].face_angle[grid.no_phi] = cos_alpha
        vertex[0].neighbour[1].face_angle[2] = cos_alpha               
        for k in range(1, grid.no_phi-1):
            temp = 1+k
            vertex[0].neighbour[temp] = FmmNeighbour()
            vertex[0].neighbour[temp].face = [k, 2+k]
            vertex[0].neighbour[temp].face_angle[k] = cos_alpha
            vertex[0].neighbour[temp].face_angle[2+k] = cos_alpha 
        vertex[0].neighbour[grid.no_phi] = FmmNeighbour()
        vertex[0].neighbour[grid.no_phi].face = [grid.no_phi-1, 1]
        vertex[0].neighbour[grid.no_phi].face_angle[grid.no_phi-1] = cos_alpha
        vertex[0].neighbour[grid.no_phi].face_angle[1] = cos_alpha
            
    elif i == grid.no_vertices-1: 
        
        # South pole
        vertex[grid.no_vertices-1].neighbour[k0] = FmmNeighbour()
        vertex[grid.no_vertices-1].neighbour[k0].face = [grid.no_vertices-2, k0+1]
        vertex[grid.no_vertices-1].neighbour[k0].face_angle[grid.no_vertices-2] = cos_alpha
        vertex[grid.no_vertices-1].neighbour[k0].face_angle[k0+1] = cos_alpha
        for k in range(1, grid.no_phi-1):
            temp = k0+k
            vertex[grid.no_vertices-1].neighbour[temp] = FmmNeighbour()
            vertex[grid.no_vertices-1].neighbour[temp].face = [temp-1, temp+1]
            vertex[grid.no_vertices-1].neighbour[temp].face_angle[temp-1] = cos_alpha
            vertex[grid.no_vertices-1].neighbour[temp].face_angle[temp+1] = cos_alpha  
        vertex[grid.no_vertices-1].neighbour[grid.no_vertices-2] = FmmNeighbour()
        vertex[grid.no_vertices-1].neighbour[grid.no_vertices-2].face = [grid.no_vertices-3, k0]
        vertex[grid.no_vertices-1].neighbour[grid.no_vertices-2].face_angle[grid.no_vertices-3] = cos_alpha
        vertex[grid.no_vertices-1].neighbour[grid.no_vertices-2].face_angle[k0] = cos_alpha
        
    else:     
        # Not a pole

        # Up neighbour (minus theta)
        if th == 1:
            j_up = 0
        else:
            j_up = grid.get_vertex_index(th-1, ph)
        vertex[i].neighbour[j_up] = FmmNeighbour()    
        
        # Down neighbour (plus theta)
        if th == grid.no_theta-2:
            j_dn = grid.no_vertices - 1
        else:
            j_dn = grid.get_vertex_index(th+1, ph)

        vertex[i].neighbour[j_dn] = FmmNeighbour() 
        
        # Left neighbour (minus phi)
        if ph == 0:
            j_lt = grid.get_vertex_index(th, grid.no_phi-1)
        else:
            j_lt = grid.get_vertex_index(th, ph-1)
        vertex[i].neighbour[j_lt] = FmmNeighbour()
        
        # Right neighbour (plus phi)
        if ph == grid.no_phi-1:
            j_rt = grid.get_vertex_index(th, 0)
        else:
            j_rt = grid.get_vertex_index(th, ph+1)
        vertex[i].neighbour[j_rt] = FmmNeighbour()
        
        # Diagonal neighbours
        if is_connect_8 == True:
            if i > grid.no_phi:
                # Up-left neighbour (minus theta, minus phi)
                if th == 1:
                    j_uplt = 0
                else:
                    if ph == 0:
                        j_uplt = grid.get_vertex_index(th-1, grid.no_phi-1)
                    else:
                        j_uplt = grid.get_vertex_index(th-1, ph-1)
                vertex[i].neighbour[j_uplt] = FmmNeighbour()
                
                # Up-right neighbour (minus theta, plus phi)
                if th == 1:
                    j_uprt = 0
                else:
                    if ph == grid.no_phi-1:
                        j_uprt = grid.get_vertex_index(th-1, 0)
                    else:
                        j_uprt = grid.get_vertex_index(th-1, ph+1)
                vertex[i].neighbour[j_uprt] = FmmNeighbour()
            
            if i < k0:
                # Down-left neighbour (plus theta, minus phi)
                if th == grid.no_theta-2:
                    j_dnlt = grid.no_vertices - 1
                else:
                    if ph == 0:
                        j_dnlt = grid.get_vertex_index(th+1, grid.no_phi-1)
                    else:
                        j_dnlt = grid.get_vertex_index(th+1, ph-1)
                vertex[i].neighbour[j_dnlt] = FmmNeighbour()
                
                # Down-right neighbour (plus theta, plus phi)
                if th == grid.no_theta-2:
                    j_dnrt = grid.no_vertices - 1
                else:
                    if ph == grid.no_phi-1:
                        j_dnrt = grid.get_vertex_index(th+1, 0)
                    else:
                        j_dnrt = grid.get_vertex_index(th+1, ph+1)
                vertex[i].neighbour[j_dnrt] = FmmNeighbour()
            
            # Find faces for 8-connectivity case
            if i <= grid.no_phi:
                # North pole adjacent
                vertex[i].neighbour[j_up].face = [j_lt, j_rt]
                vertex[i].neighbour[j_up].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_up].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_rt].face = [j_up, j_dnrt]
                vertex[i].neighbour[j_rt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_rt].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face = [j_rt, j_dn]
                vertex[i].neighbour[j_dnrt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dn].face = [j_dnlt, j_dnrt]
                vertex[i].neighbour[j_dn].face_angle[j_dnlt] = cos_alpha
                vertex[i].neighbour[j_dn].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnlt].face = [j_dn, j_lt]
                vertex[i].neighbour[j_dnlt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dnlt].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_lt].face = [j_up, j_dnlt]
                vertex[i].neighbour[j_lt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_lt].face_angle[j_dnlt] = cos_alpha
                
            elif i >= k0:
                # South pole adjacent
                vertex[i].neighbour[j_up].face = [j_uplt, j_uprt]
                vertex[i].neighbour[j_up].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_up].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_uprt].face = [j_up, j_rt]
                vertex[i].neighbour[j_uprt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uprt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_rt].face = [j_uprt, j_dn]
                vertex[i].neighbour[j_rt].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_rt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dn].face = [j_lt, j_rt]
                vertex[i].neighbour[j_dn].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_dn].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_lt].face = [j_dn, j_uplt]
                vertex[i].neighbour[j_lt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_lt].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_uplt].face = [j_up, j_lt]
                vertex[i].neighbour[j_uplt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uplt].face_angle[j_lt] = cos_alpha
                
            else:
                # Not pole adjacent
                vertex[i].neighbour[j_up].face = [j_uplt, j_uprt]
                vertex[i].neighbour[j_up].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_up].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_uprt].face = [j_up, j_rt]
                vertex[i].neighbour[j_uprt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uprt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_rt].face = [j_uprt, j_dnrt]
                vertex[i].neighbour[j_rt].face_angle[j_uprt] = cos_alpha
                vertex[i].neighbour[j_rt].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face = [j_rt, j_dn]
                vertex[i].neighbour[j_dnrt].face_angle[j_rt] = cos_alpha
                vertex[i].neighbour[j_dnrt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dn].face = [j_dnlt, j_dnrt]
                vertex[i].neighbour[j_dn].face_angle[j_dnlt] = cos_alpha
                vertex[i].neighbour[j_dn].face_angle[j_dnrt] = cos_alpha
                vertex[i].neighbour[j_dnlt].face = [j_dn, j_lt]
                vertex[i].neighbour[j_dnlt].face_angle[j_dn] = cos_alpha
                vertex[i].neighbour[j_dnlt].face_angle[j_lt] = cos_alpha
                vertex[i].neighbour[j_lt].face = [j_uplt, j_dnlt]
                vertex[i].neighbour[j_lt].face_angle[j_uplt] = cos_alpha
                vertex[i].neighbour[j_lt].face_angle[j_dnlt] = cos_alpha
                vertex[i].neighbour[j_uplt].face = [j_up, j_lt]
                vertex[i].neighbour[j_uplt].face_angle[j_up] = cos_alpha
                vertex[i].neighbour[j_uplt].face_angle[j_lt] = cos_alpha
            
        else: 
            # Find faces for 4-connectivity case
            vertex[i].neighbour[j_up].face = [j_lt, j_rt]
            vertex[i].neighbour[j_up].face_angle[j_lt] = cos_alpha
            vertex[i].neighbour[j_up].face_angle[j_rt] = cos_alpha
            vertex[i].neighbour[j_dn].face = [j_lt, j_rt]
            vertex[i].neighbour[j_dn].face_angle[j_lt] = cos_alpha
            vertex[i].neighbour[j_dn].face_angle[j_rt] = cos_alpha
            vertex[i].neighbour[j_lt].face = [j_up, j_dn]
            vertex[i].neighbour[j_lt].face_angle[j_up] = cos_alpha
            vertex[i].neighbour[j_lt].face_angle[j_dn] = cos_alpha
            vertex[i].neighbour[j_rt].face = [j_up, j_dn]
            vertex[i].neighbour[j_rt].face_angle[j_up] = cos_alpha
            vertex[i].neighbour[j_rt].face_angle[j_dn] = cos_alpha

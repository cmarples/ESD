"""! 
@brief Defines classes used in the fast marching method.
@file classes.py
@author Callum Marples
- Created by Callum Marples on 12/12/2022.
- Last modified on 02/02/2022.
"""

from math import sqrt, inf
import numpy as np

class FmmVertex:
    """! The vertex class.
    Contains information for a given vertex in a mesh used for the
    Dijkstra's algorithm or the fast marching method.
    """
    def __init__(self, index, carts=np.zeros(3), is_endpoint=False):
        """! The FmmVertex initialiser.
        @param index : int \n
            The scalar index of the vertex.
        @param carts : 3-element NumPy array (optional) \n
            Cartesian coordinates of the vertex in 3 dimensional space. Defaults to [0.0, 0.0, 0.0].
        @param is_endpoint : bool (optional) \n
            If True, this vertex represents an end point. Defaults to False.
        """
        ## Scalar index of the vertex.
        self.index = index
        ## Cartesian coordinates of the vertex.
        self.carts = carts
        ## Dictionary containing all neighbours of the vertex
        self.neighbour = {}
        ## Set of all neighbours of the vertex
        self.neighbour_set = set()
        ## Integer giving the icosahedral face the vertex belongs to (used for building an icosahedral mesh)
        self.ico = []

class FmmNeighbour:
    """! The neighbour class.
    Contains information for a given neighbour, j, of vertex i. This information 
    includes Euclidean distance between the neighbours, the faces to which the 
    pair belong and the angles in these faces.
    """
    def __init__(self):
        """! The FmmNeighbour initialiser.
        """
        ## Euclidean distance between neighbours.
        self.distance = -1.0
        ## List of two integers, giving the common neighbours of the pair stored here.
        ## This gives the two faces to which this edge belongs.
        self.face = []
        ## Dictionary giving the angles at the host vertex, in the two faces specified by face.
        self.face_angle = {}

class FmmMesh:
    """! The mesh class.
    Contains a list of FmmVertex objects along with any other information
    pertaining to the mesh.
    """
    def __init__(self):
        """! The FmmMesh initialiser.
        """
        ## List of vertices in the mesh.
        self.vertex = []
        ## Number of vertices in the mesh.
        self.no_vertices = 0
        ## Polar ("pol") or Icosahedral ("ico").
        self.type = "none"
        
    def precalculate(self):
        """! Precalculates all edge lengths and face angles on a mesh.
        @see gen_pol_mesh
        @see gen_ico_mesh
        """
        ## Smallest edge length.
        self.min_edge = inf
        ## Largest edge length.
        self.max_edge = 0.0
        for i in range(self.no_vertices):
            for j in self.vertex[i].neighbour:
                if i > j:
                    self.vertex[i].neighbour[j].distance = sqrt( (self.vertex[i].carts[0]-self.vertex[j].carts[0])**2.0 +
                                                            (self.vertex[i].carts[1]-self.vertex[j].carts[1])**2.0 +
                                                            (self.vertex[i].carts[2]-self.vertex[j].carts[2])**2.0 )
                    self.vertex[j].neighbour[i].distance = self.vertex[i].neighbour[j].distance
                    # Check if this distance is a new maximum or minimum.
                    if self.vertex[i].neighbour[j].distance > self.max_edge:
                        self.max_edge = self.vertex[i].neighbour[j].distance
                    if self.vertex[i].neighbour[j].distance < self.min_edge:
                        self.min_edge = self.vertex[i].neighbour[j].distance
                        
        ## Number of obtuse angles in the mesh.
        self.no_obtuse = 0
        ## Cosine of the smallest face angle.
        self.min_angle = -2.0
        ## Cosine of the largest face angle.
        self.max_angle = 2.0
        for i in range(self.no_vertices):        
            for j in self.vertex[i].neighbour:
                for k in self.vertex[i].neighbour[j].face_angle.keys():           
                    if self.vertex[i].neighbour[j].face_angle[k] == 2.0:
                        w1 = self.vertex[j].carts - self.vertex[i].carts
                        w2 = self.vertex[k].carts - self.vertex[i].carts
                        a = self.vertex[i].neighbour[j].distance
                        b = self.vertex[i].neighbour[k].distance
                        self.vertex[i].neighbour[j].face_angle[k] = (w1[0]*w2[0] + w1[1]*w2[1] + w1[2]*w2[2]) / (a*b)
                        self.vertex[i].neighbour[k].face_angle[j] = self.vertex[i].neighbour[j].face_angle[k]
                        # Check if this angle is obtuse, a new maximum or a new minimum.
                        if self.vertex[i].neighbour[j].face_angle[k] < 0.0: # Obtuse angle (cosine negative)
                            self.no_obtuse += 1
                        if self.vertex[i].neighbour[j].face_angle[k] < self.max_angle:
                            self.max_angle = self.vertex[i].neighbour[j].face_angle[k]
                        if self.vertex[i].neighbour[j].face_angle[k] > self.min_angle:
                            self.min_angle = self.vertex[i].neighbour[j].face_angle[k]
                        
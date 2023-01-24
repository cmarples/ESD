# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 20:08:13 2022

@author: Cal
"""

import math
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
        self.index = index
        self.is_endpoint = is_endpoint
        self.carts = carts
        self.neighbour = {}
        self.neighbour_set = set()
        self.ico = []

class FmmNeighbour:
    """! The neighbour class.
    Contains information for a given neighbour, j, of vertex i.
    """
    def __init__(self):
        """! The FmmNeighbour initialiser.
        """
        self.distance = -1.0
        self.face = []
        self.face_angle = {}
        self.ico_face = []

class FmmMesh:
    """! The mesh class.
    Contains a list of FmmVertex objects along with any other information
    pertaining to the mesh.
    """
    def __init__(self):
        """! The FmmMesh initialiser.
        """
        self.vertex = []
        self.no_vertices = 0
        self.type = "none"
        
    def precalculate(self):
        """! Precalculates all edge lengths and face angles on a mesh.
        @see gen_pol_mesh
        @see gen_ico_mesh
        """
        for i in range(self.no_vertices):
            for j in self.vertex[i].neighbour:
                if i > j:
                    self.vertex[i].neighbour[j].distance = math.sqrt( (self.vertex[i].carts[0]-self.vertex[j].carts[0])**2.0 +
                                                            (self.vertex[i].carts[1]-self.vertex[j].carts[1])**2.0 +
                                                            (self.vertex[i].carts[2]-self.vertex[j].carts[2])**2.0 )
                    self.vertex[j].neighbour[i].distance = self.vertex[i].neighbour[j].distance
                    
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
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 20:08:13 2022

@author: Cal
"""

import math
import numpy as np

class FmmVertex:
    """! The vertex class.
    Contains information for a given vertex in a grid used for the
    Dijkstra's algorithm or the fast marching method.
    """
    def __init__(self, index, carts=np.zeros(3), is_endpoint=False):
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
        self.distance = -1.0
        self.face = []
        self.face_angle = {}
        self.ico_face = []

class FmmGrid:
    """! The grid class.
    Contains a list of FmmVertex objects along with any other information
    pertaining to the grid.
    """
    def __init__(self):
        self.vertex = []
        self.no_vertices = 0
        self.type = "none"
        
    def precalculate(self):

        for i in range(self.no_vertices):
            for j in self.vertex[i].neighbour:
                if i > j:
                    # Calculate distance between vertices i and j
                    #dist = get_distance(vertex, i, j)
    
                    self.vertex[i].neighbour[j].distance = math.sqrt( (self.vertex[i].carts[0]-self.vertex[j].carts[0])**2.0 +
                                                            (self.vertex[i].carts[1]-self.vertex[j].carts[1])**2.0 +
                                                            (self.vertex[i].carts[2]-self.vertex[j].carts[2])**2.0 )
                    self.vertex[j].neighbour[i].distance = self.vertex[i].neighbour[j].distance
                    
        for i in range(self.no_vertices):        
            for j in self.vertex[i].neighbour:
                for k in self.vertex[i].neighbour[j].face_angle.keys():
                    #cos_alpha = get_angle(vertex, i, j, k)
                    
                    if self.vertex[i].neighbour[j].face_angle[k] == 2.0:
                        w1 = self.vertex[j].carts - self.vertex[i].carts
                        w2 = self.vertex[k].carts - self.vertex[i].carts
                        a = self.vertex[i].neighbour[j].distance
                        b = self.vertex[i].neighbour[k].distance
                        self.vertex[i].neighbour[j].face_angle[k] = (w1[0]*w2[0] + w1[1]*w2[1] + w1[2]*w2[2]) / (a*b)
                        self.vertex[i].neighbour[k].face_angle[j] = self.vertex[i].neighbour[j].face_angle[k]
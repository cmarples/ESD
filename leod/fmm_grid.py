# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:36:12 2022

@author: Callum Marples

Grid class for use in the Fast Marching Method.

"""

# This class contains the set of all FmmVertex objects for a given surface.
# It also includes the neighbours and their distances
class FmmGrid:
    # Constructor
    def __init__(self, is_Dijkstra=False):
        self.is_Dijkstra = is_Dijkstra
        self.vertex = []

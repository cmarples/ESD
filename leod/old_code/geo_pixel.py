# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:06:00 2022

@author: Callum Marples
"""

import numpy as np

# This class contains information for a pixel in the grid used in the fast
# marching method. An array of these objects is associated with a particular
# EllipsoidShape object.
class GeoPixel:
    def __init__(self, pixel_index, theta_index, phi_index, no_pixels, carts=np.zeros(3), is_endpoint=False, is_refine=False):
        self.pixel_index = pixel_index
        self.theta_index = theta_index
        self.phi_index = phi_index
        self.set_pole(no_pixels)
        self.is_endpoint = is_endpoint
        if is_refine == False:
            self.carts = np.array(carts)
            self.neighbour = []
            self.neighbour_distance = []
            self.is_border = False
            

    def set_pole(self, no_pixels):
        if self.pixel_index == 0:
            self.phi_index = 0
            self.is_north = True
            self.is_south = False
        elif self.pixel_index == no_pixels-1:
            self.phi_index = 0
            self.is_north = False
            self.is_south = True
        else:
            self.is_north = False
            self.is_south = False
    
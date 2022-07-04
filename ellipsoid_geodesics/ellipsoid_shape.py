# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 16:13:12 2022

@author: Callum Marples
"""

class EllipsoidShape:
    
    def __init__(self, a=1.0, b=1.0, c=1.0):
        self.set_axes(a, b, c)

    def is_sphere(self):
        if self.a_axis == self.b_axis and self.b_axis == self.c_axis:
            return True
        else:
            return False
        
    def set_axes(self, a, b, c):
        # Check that the inputs are valid
        if not (isinstance(a, (int, float)) and a > 0):
             raise ValueError(f"positive a-axis expected, got {a}")
        if not (isinstance(b, (int, float)) and b > 0):
             raise ValueError(f"positive b-axis expected, got {b}")
        if not (isinstance(c, (int, float)) and c > 0):
             raise ValueError(f"positive c-axis expected, got {c}")
        # Set axis length values
        self.a_axis = a
        self.b_axis = b
        self.c_axis = c
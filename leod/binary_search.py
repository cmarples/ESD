# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 13:42:23 2022

@author: Callum Marples
"""

import math

# Find index of closest value in list v to number x,
# using the binary search algorithm.
# It is assumed that the elements in v are monotonically increasing.
def binary_search(v, x):
    # Initialise
    L = 0
    R = len(v)
    # Find index L such that : v[L] <= x < v[L+1]
    while L != R:
        M = math.ceil( 0.5*(L + R) )
        if x < v[M]:
            R = M - 1
        else:
            L = M
    # Output the index whose value is closest to x
    dL = x - v[L]
    dR = v[L+1] - x
    if dL < dR:
        return L
    else:
        return L + 1

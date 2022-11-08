# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 14:59:44 2022

@author: Callum Marples

Given ellipsoid (a,b,c) and a sphere of radius r (centred at an ellipsoid
surface point), find a sample of points on the intersecting curve.

Using the method of Yu et al.
"""

import math
from scipy import optimize

from leod.newton_raphson import newton_raphson

# The function that implicit defines of the curve Cp
def h_implicit(th, ph, centre, a, b, c, r2):
    s_th = math.sin(th)
    return ( (a*s_th*math.cos(ph) - centre[0])**2.0 + 
             (b*s_th*math.sin(ph) - centre[1])**2.0 +
             (c*math.cos(th) - centre[2])**2.0 - r2 )
    
# Obtain an initial point on parametric space curve Cp
def initial_point(centre, r, alpha, beta, c, u0):
    # Define function and its derivative
    h = lambda u: ( (alpha*math.sqrt(1-u*u) - centre[0])**2.0 +
                    (beta*math.sqrt(1-u*u) - centre[1])**2.0 +
                    (c*u - centre[2])**2.0 - r*r )
    h_prime = lambda u: 2.0 * ( c*(c*u - centre[2])  - 
                                u*(alpha*(alpha*math.sqrt(1-u*u) - centre[0]) + beta*(beta*math.sqrt(1-u*u) - centre[1])) / math.sqrt(1-u*u) )
    #u = newton_raphson(h, h_prime, u0)
    sol = optimize.root_scalar(h, bracket=[0.0, 1.0], method='brentq')
    
    return math.acos(sol.root)


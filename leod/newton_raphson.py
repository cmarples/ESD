# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 15:05:32 2022

@author: Callum Marples

Newton-Raphson method for finding roots of functions of one variable.
"""

from math import fabs

# Evaluate root of a function using the Newton-Raphson method.
# No safeguarding - applicable when root bounds are unknown or infinite.
def newton_raphson(fun, fun_prime, x, tol=1e-14):
    x_diff = 1.0
    f_diff = 1.0
    f = fun(x)
    while fabs(x_diff) > tol and fabs(f_diff) > tol:
        f_prime = fun_prime(x)
        x_diff = f / f_prime
        x -= x_diff
        f_new = fun(x)
        f_diff = f_new - f
        f = f_new
    return x
        
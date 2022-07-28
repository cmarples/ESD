# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 19:45:31 2022

@author: Callum Marples

Unit sphere

Calculate example distance for x_ref = 1 over a range of no_div and plot

Is there a sudden jump at no_div = 13?
"""

import os
os.chdir("..")

import math
import numpy as np
import csv
import time

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid
from leod.geo_fmm import GeoFMM

E = EllipsoidShape(1.0, 1.0, 1.0)

th_0 = 90.0 * math.pi / 180.0
ph_0 = 0.0
th_1 = 50.0 * math.pi / 180.0
ph_1 = 60.0 * math.pi / 180.0

no_div = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
d = [0.0] * 12
G = GeoGrid(E, 200, 200)
for i in range(12):
    tic = time.perf_counter()
    F = GeoFMM(G, th_0, ph_0)
    d[i] = F.calculate_geodesics(2, th_1, ph_1, is_refine=True, refine_range=1, 
                                 refine_theta=no_div[i], refine_phi=no_div[i])
    del F
    toc = time.perf_counter()
    print(toc - tic)

file_name = 'dist_vs_divisions.txt'
with open(file_name, mode="w", newline='') as f:
    
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(d)    
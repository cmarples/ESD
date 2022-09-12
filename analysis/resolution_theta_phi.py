# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 15:40:06 2022

@author: Callum Marples
"""

import math
import numpy as np
import time

from leod.ellipsoid_shape import EllipsoidShape
from leod.geo_grid import GeoGrid
from leod.geo_fmm import GeoFMM

E = EllipsoidShape(1.0, 1.0, 1.0)
th_0 = 90.0 * math.pi / 180.0
ph_0 = 0.0
th_1 = 50.0 * math.pi / 180.0
ph_1 = 60.0 * math.pi / 180.0

n_res = np.arange(10, 370, step=10)
m = len(n_res)

x_ref = 10
n_div = 5

distances = np.zeros([m, m])
dist_rfnd = np.zeros([m, m])
times = np.zeros([m, m])
times_rfnd = np.zeros([m, m])


for i in range(m):        # i : theta index
    for j in range(m):    # j : phi index
        
        # Unrefined
        G1 = GeoGrid(E, n_res[i], n_res[j])
        F1 = GeoFMM(G1, th_0, ph_0)
        tic = time.perf_counter()
        distances[i][j] = F1.calculate_geodesics(2, th_1, ph_1, False)
        toc = time.perf_counter()
        times[i][j] = toc - tic
        
        # Refined
        G2 = GeoGrid(E, n_res[i], n_res[j])
        F2 = GeoFMM(G2, th_0, ph_0)
        tic = time.perf_counter()
        dist_rfnd[i][j] = F2.calculate_geodesics(2, th_1, ph_1, True, x_ref, n_div, n_div)
        toc = time.perf_counter()
        times_rfnd[i][j] = toc - tic
        
    if (i+1)%6 == 0:
        print(['i = ', i, ' out of 35 complete.'])


# Write to files
np.savetxt('data/res_theta_phi.csv', distances)
np.savetxt('data/res_theta_phi_refined.csv', dist_rfnd)
np.savetxt('data/res_theta_phi_times.csv', times)
np.savetxt('data/res_theta_phi_refined_times.csv', times_rfnd)

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 21:51:18 2022

@author: Callum Marples
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import os


    
dist = np.loadtxt('../data/res_theta_phi.csv', delimiter = ",")
time = np.loadtxt('../data/res_theta_phi_times.csv', delimiter = ",")
dist_rfnd = np.loadtxt('../data/res_theta_phi_refined.csv', delimiter = ",")
time_rfnd = np.loadtxt('../data/res_theta_phi_refined_times.csv', delimiter = ",")

os.chdir("../..")
from leod.sphere_geodesics import great_circle_distance
# True distance
th_0 = 90.0 * math.pi / 180.0
ph_0 = 0.0

th_1 = 50.0 * math.pi / 180.0
ph_1 = 60.0 * math.pi / 180.0
    
s = great_circle_distance(1.0, th_0, ph_0, th_1, ph_1)

d1 = dist - s

# Subset with low resolution removed
n_res = np.arange(70, 370, step=10)
m = len(n_res)
d = np.zeros([m, m])
d2 = np.zeros([m, m])
for i in range(m):        # i : theta index
    for j in range(m):    # j : phi index
        d[i][j] = dist[i+6][j+6]
        d2[i][j] = dist_rfnd[i+6][j+6]

fig1, ax = plt.subplots(1, 1, figsize=(10, 10))
# Plot difference image in 2d
map_name = "cool"
plt.imshow(dist - s, cmap=plt.get_cmap(map_name), origin='lower') 

xt = np.arange(60, 420, 60)
xx = np.arange(5, 36, 6)
plt.xticks(xx, xt)
yt = np.arange(60, 420, 60)
yy = np.arange(5, 36, 6)
plt.yticks(yy, yt)

plt.xlabel('$n_\phi$')
plt.ylabel('$n_\\theta$')
    
plt.colorbar()

fig2, ax = plt.subplots(1, 1, figsize=(10, 10))
# Plot difference image in 2d
map_name = "cool"
plt.imshow(dist_rfnd - s, cmap=plt.get_cmap(map_name), origin='lower') 

xt = np.arange(60, 420, 60)
xx = np.arange(5, 36, 6)
plt.xticks(xx, xt)
yt = np.arange(60, 420, 60)
yy = np.arange(5, 36, 6)
plt.yticks(yy, yt)

plt.xlabel('$n_\phi$')
plt.ylabel('$n_\\theta$')
    
plt.colorbar()


fig3, ax = plt.subplots(1, 1, figsize=(10, 10))
plt.imshow(d-s, cmap=plt.get_cmap(map_name), origin='lower') 

xt = np.arange(120, 420, 60)
xx = np.arange(5, 30, 6)
plt.xticks(xx, xt)
yt = np.arange(120, 420, 60)
yy = np.arange(5, 30, 6)
plt.yticks(yy, yt)

plt.xlabel('$n_\phi$')
plt.ylabel('$n_\\theta$')
    
plt.colorbar()


fig4, ax = plt.subplots(1, 1, figsize=(10, 10))
plt.imshow(d2-s, cmap=plt.get_cmap(map_name), origin='lower') 

xt = np.arange(120, 420, 60)
xx = np.arange(5, 30, 6)
plt.xticks(xx, xt)
yt = np.arange(120, 420, 60)
yy = np.arange(5, 30, 6)
plt.yticks(yy, yt)

plt.xlabel('$n_\phi$')
plt.ylabel('$n_\\theta$')
    
plt.colorbar()
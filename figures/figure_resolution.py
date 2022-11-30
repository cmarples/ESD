# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 17:57:40 2022

@author: Callum Marples

Plot fast marching resolution data.
Data generated in:
    - apps/run_geodesic_resolution.py
    - apps/run_geodesic_resolution_comparison.py
Note that array sizes are hardcoded in this script.
"""

import math
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
from figure_size import set_size

os.chdir("..")

save_on = True
plt.style.use('tex_symmetry')
fig_width = 345.0

# Read shape and alternate method data.
a = [0] * 5
b = [0] * 5
c = [0] * 5
d_true = ['N/A'] * 5
d_geo = ['N/A'] * 5
d_bvm = ['N/A'] * 5
d_taxi = ['N/A'] * 5
rows = ['N/A'] * 5
file_comparison = 'data/resolution/resolution_comparison.txt'
with open(file_comparison) as f:
    csv_reader = csv.reader(f, delimiter=',')
    i = 0
    for row in csv_reader:
        if i > 0 and i < 6:
            rows[i-1] = row
        i += 1
for i in range(5):
    if rows[i][4] != "N/A":
        d_true[i] = float(rows[i][4])
    if rows[i][5] != "N/A":
        d_geo[i] = float(rows[i][5])
    if rows[i][6] != "N/A":
        d_bvm[i] = float(rows[i][6])
    if rows[i][7] != "N/A":
        d_taxi[i] = float(rows[i][7])

# Read number of vertices data
no_vertices = [0] * 25
no_vertices_polar = [0] * 25
no_vertices_icosahedral = [0] * 25
shape = [0] * 5
file_number = 'data/resolution/resolution_number_of_vertices.txt'
with open(file_number) as f:
    csv_reader = csv.reader(f, delimiter=',')
    i = 0
    for row in csv_reader:
        if i < 25:
            no_vertices[i] = int(row[0])
            no_vertices_polar[i] = int(row[1])
            no_vertices_icosahedral[i] = int(row[2])
        i += 1

# Read resolution data
d_dijkstra_4 = np.zeros([len(no_vertices), len(shape)])
d_fmm_4 = np.zeros([len(no_vertices), len(shape)])
d_dijkstra_8 = np.zeros([len(no_vertices), len(shape)])
d_fmm_8 = np.zeros([len(no_vertices), len(shape)])
d_dijkstra_tri = np.zeros([len(no_vertices), len(shape)])
d_fmm_tri = np.zeros([len(no_vertices), len(shape)])
d_dijkstra_split = np.zeros([len(no_vertices), len(shape)])
d_fmm_split = np.zeros([len(no_vertices), len(shape)])


t_polar_gen_4 = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_4 = np.zeros([len(no_vertices), len(shape)])
t_fmm_4 = np.zeros([len(no_vertices), len(shape)])


t_polar_gen_8 = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_8 = np.zeros([len(no_vertices), len(shape)])
t_fmm_8 = np.zeros([len(no_vertices), len(shape)])


t_tri_gen = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_tri = np.zeros([len(no_vertices), len(shape)])
t_fmm_tri = np.zeros([len(no_vertices), len(shape)])


t_tri_split_gen = np.zeros([len(no_vertices), len(shape)])
t_dijkstra_split = np.zeros([len(no_vertices), len(shape)])
t_fmm_split = np.zeros([len(no_vertices), len(shape)])


shape_name = ['sphere', 'wgs84', 'oblate', 'prolate', 'triaxial']
for j in range(len(shape)):
    file_shape = 'data/resolution/resolution_' + shape_name[j] + '.txt'
    rows = ['N/A'] * 25
    with open(file_shape) as f:
        csv_reader = csv.reader(f, delimiter=',')
        i = 0
        for row in csv_reader:
            if i > 0 and i < 26:
                rows[i-1] = row
            i += 1
    # Extract data and convert to floating point
    for i in range(len(no_vertices)):
        d_dijkstra_4[i][j] = rows[i][0]
        d_fmm_4[i][j] = rows[i][1]
        d_dijkstra_8[i][j] = rows[i][2]
        d_fmm_8[i][j] = rows[i][3]
        d_dijkstra_tri[i][j] = rows[i][4]
        d_fmm_tri[i][j] = rows[i][5]
        d_dijkstra_split[i][j] = rows[i][6]
        d_fmm_split[i][j] = rows[i][7]


fig = [0] * 5
ax = [0] * 5

# Plot distance data
for j in range(5):
    
    fig[j], ax[j] = plt.subplots(1, 1, figsize=set_size(fig_width, height=1.0))  
    
    ax[j].scatter(np.log10(no_vertices_polar), d_dijkstra_4[:, j], c="none", edgecolor="c", marker="s")
    ax[j].scatter(np.log10(no_vertices_polar), d_dijkstra_8[:, j], c="none", edgecolor="b", marker="o")
    
    
    ax[j].scatter(np.log10(no_vertices_polar), d_fmm_4[:, j], c="c", edgecolor="c", marker="s")
    ax[j].scatter(np.log10(no_vertices_polar), d_fmm_8[:, j], c="b", edgecolor="b", marker="o")
    
    if j < 2: # Sphere or WGS84
        ax[j].scatter(np.log10(no_vertices_polar), d_dijkstra_tri[:, j], c="none", edgecolor="k", marker="^")
        ax[j].scatter(np.log10(no_vertices_polar), d_fmm_tri[:, j], c="k", edgecolor="k", marker="^")
    else:
        # Plot icosahedral triangulation distances, avoiding repetition when splitting has no effect.
        for i in range(25):
            ax[j].scatter(math.log10(no_vertices_polar[i]), d_dijkstra_split[i, j], c="none", edgecolor="k", marker="v")
            ax[j].scatter(np.log10(no_vertices_polar[i]), d_fmm_split[i, j], c="k", edgecolor="k", marker="v")
            if d_dijkstra_tri[i, j] != d_dijkstra_split[i, j]:
                ax[j].scatter(math.log10(no_vertices_polar[i]), d_dijkstra_split[i, j], c="none", edgecolor="k", marker="v")
            if d_fmm_tri[i, j] != d_fmm_split[i, j]:    
                ax[j].scatter(np.log10(no_vertices_polar[i]), d_fmm_tri[i, j], c="k", edgecolor="k", marker="^")
    
    # Plot alternate distances.            
    if j == 0:
        ax[j].plot(np.log10(no_vertices), [d_true[j]]*25, 'k', linestyle="dotted")
    elif j > 0:
        if j < 4:
            ax[j].plot(np.log10(no_vertices), [d_geo[j]]*25, 'k', linestyle="solid")
        ax[j].plot(np.log10(no_vertices), [d_bvm[j]]*25, 'k', linestyle="dashed")
    if j < 4:
        ax[j].plot(np.log10(no_vertices), [d_taxi[j]]*25, 'k', linestyle="dashdot")
    # Set x-axis tick labels (number of vertices).  
    ax[j].set_xticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    ax[j].set_xticklabels(["$10^2$", "$10^{2.5}$", "$10^3$", "$10^{3.5}$", "$10^4$", "$10^{4.5}$", "$10^5$"])    
        
       
        
       
        
        
    
    
    






      
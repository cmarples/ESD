# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 13:04:26 2022

@author: Callum Marples

Plot fast marching time vs. resolution data.
Data generated in:
    - apps/run_geodesic_resolution.py
    - apps/run_geodesic_resolution_comparison.py
Note that array sizes are hardcoded in this script.
"""

import numpy as np
import csv
import os
import matplotlib.pyplot as plt
from figure_size import set_size

os.chdir("..")

save_on = True
plt.style.use('tex')
fig_width = 345.0


# Read shape and alternate method data.
a = [0] * 5
b = [0] * 5
c = [0] * 5

t_true = ['N/A'] * 5
t_geo = ['N/A'] * 5
t_bvm = ['N/A'] * 5
t_taxi = ['N/A'] * 5
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
        t_true[i] = float(rows[i][8])
    if rows[i][5] != "N/A":
        t_geo[i] = float(rows[i][9])
    if rows[i][6] != "N/A":
        t_bvm[i] = float(rows[i][10])
    if rows[i][7] != "N/A":
        t_taxi[i] = float(rows[i][11])
        


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
        t_polar_gen_4[i][j] = rows[i][8]
        t_dijkstra_4[i][j] = rows[i][9]
        t_fmm_4[i][j] = rows[i][10]
        t_polar_gen_8[i][j] = rows[i][11]
        t_dijkstra_8[i][j] = rows[i][12]
        t_fmm_8[i][j] = rows[i][13]
        t_tri_gen[i][j] = rows[i][14]
        t_dijkstra_tri[i][j] = rows[i][15]
        t_fmm_tri[i][j] = rows[i][16]
        t_tri_split_gen[i][j] = rows[i][17]
        t_dijkstra_split[i][j] = rows[i][18]
        t_fmm_split[i][j] = rows[i][19]

fig_dist = [0] * 6
ax_dist = [0] * 6

m = [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]
ms = 9
lw = 1
og = [[1, 0.65, 0]]

fig_gen, ax_gen= plt.subplots(3, 2, figsize=set_size(fig_width, height=1.0, subplots=(3,2)))
fig_fmm, ax_fmm = plt.subplots(3, 2, figsize=set_size(fig_width, height=1.0, subplots=(3,2)))
shape_name = ["Sphere", "WGS84", "Oblate", "Prolate", "Triaxial"]
shape_label = ["(a)", "(b)", "(c)", "(d)", "(e)"]

# Plot distance data
for j in range(5): 
    
    # Plot non-wavefront distances.            
    if j == 0:
        ax_gen[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [t_true[j]]*25, 'k', linestyle="dotted", linewidth=1.5)
    elif j > 0:
        if j < 4:
            ax_gen[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [t_geo[j]]*25, 'k', linestyle="solid", linewidth=lw)
        ax_gen[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [t_bvm[j]]*25, 'k', linestyle="dashed", linewidth=lw)
    
    ax_gen[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_polar_gen_4[:, j], c="c", edgecolor="c", marker="s", s=ms)    
    ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_fmm_4[:, j], c="c", edgecolor="c", marker="s", s=ms)
    
    ax_gen[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_polar_gen_8[:, j], c="b", edgecolor="b", marker="o", s=ms)    
    ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_fmm_8[:, j], c="b", edgecolor="b", marker="o", s=ms)
    
    ax_gen[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_tri_gen[:, j], c=og, edgecolor=og, marker="^", s=ms)
    ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_fmm_tri[:, j], c=og, edgecolor=og, marker="^", s=ms)
    
    ax_gen[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_tri_split_gen[:, j], c="r", edgecolor="r", marker="v", s=ms)
    ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), t_fmm_split[:, j], c="r", edgecolor="r", marker="v", s=ms)        
            
    # Set x-axis tick labels (number of vertices).  
    ax_fmm[m[j][0]][m[j][1]].set_xticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    ax_fmm[m[j][0]][m[j][1]].set_xticklabels(["$10^2$", "$10^{2.5}$", "$10^3$", "$10^{3.5}$", "$10^4$", "$10^{4.5}$", "$10^5$"]) 
    ax_gen[m[j][0]][m[j][1]].set_xticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    ax_gen[m[j][0]][m[j][1]].set_xticklabels(["$10^2$", "$10^{2.5}$", "$10^3$", "$10^{3.5}$", "$10^4$", "$10^{4.5}$", "$10^5$"])    
      
    # Axis labels.
    ax_fmm[m[j][0]][m[j][1]].set_xlabel('Number of vertices')
    ax_fmm[m[j][0]][m[j][1]].set_ylabel('Run-time (s)')
    ax_gen[m[j][0]][m[j][1]].set_xlabel('Number of vertices')
    ax_gen[m[j][0]][m[j][1]].set_ylabel('Run-time (s)')
    
    # Titles and subfigures.
    ax_fmm[m[j][0]][m[j][1]].text(0.5, -0.42, shape_label[j], ha="center",  transform=ax_fmm[m[j][0]][m[j][1]].transAxes)
    ax_gen[m[j][0]][m[j][1]].text(0.5,-0.42, shape_label[j], ha="center",  transform=ax_gen[m[j][0]][m[j][1]].transAxes)
    ax_fmm[m[j][0]][m[j][1]].text(0.5, 1.05, shape_name[j], ha="center",  transform=ax_fmm[m[j][0]][m[j][1]].transAxes)
    ax_gen[m[j][0]][m[j][1]].text(0.5, 1.05, shape_name[j], ha="center",  transform=ax_gen[m[j][0]][m[j][1]].transAxes)
    
    ax_fmm[m[j][0]][m[j][1]].spines['top'].set_visible(False)  
    ax_fmm[m[j][0]][m[j][1]].spines['right'].set_visible(False) 
    ax_gen[m[j][0]][m[j][1]].spines['top'].set_visible(False)  
    ax_gen[m[j][0]][m[j][1]].spines['right'].set_visible(False) 
    
    ax_fmm[m[j][0]][m[j][1]].set_ylim(0, 2)
    ax_gen[m[j][0]][m[j][1]].set_ylim(0, 17)
    
# Add legend in final ax slot.
x = -1
y = -1
j = 5

ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c="c", edgecolor="c", marker="s", s = ms, label = "FMM (polar 4)")
ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c="b", edgecolor="b", marker="o", s = ms, label = "FMM (polar 8)")
ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c=og, edgecolor=og, marker="^", s = ms, label = "FMM (icosahedral)")
ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c="r", edgecolor="r", marker="v", s = ms, label = "FMM (split icosahedral)")

ax_gen[m[j][0]][m[j][1]].scatter(x, y, c="c", edgecolor="c", marker="s", s = ms, label = "4-neighbour polar")
ax_gen[m[j][0]][m[j][1]].scatter(x, y, c="b", edgecolor="b", marker="o", s = ms, label = "8-neighbour polar")
ax_gen[m[j][0]][m[j][1]].scatter(x, y, c=og, edgecolor=og, marker="^", s = ms, label = "Icosahedral triangulation")
ax_gen[m[j][0]][m[j][1]].scatter(x, y, c="r", edgecolor="r", marker="v", s = ms, label = "Icosahedral with splitiing")

ax_gen[m[j][0]][m[j][1]].plot([-2, -1], [-2, -1], 'k', linestyle="dotted", label = "True distance")
ax_gen[m[j][0]][m[j][1]].plot(x, y, 'k', linestyle="solid", label = "GeographicLib")
ax_gen[m[j][0]][m[j][1]].plot(x, y, 'k', linestyle="dashed", label = "Boundary value method")

  
ax_fmm[m[j][0]][m[j][1]].legend(loc="center left")
ax_fmm[m[j][0]][m[j][1]].axis('off')  
ax_gen[m[j][0]][m[j][1]].legend(loc="center left")
ax_gen[m[j][0]][m[j][1]].axis('off') 
 
plt.figure(fig_gen.number)  
plt.legend(frameon=False)
plt.xlim([0, 1])  
plt.ylim([0, 1]) 
plt.subplots_adjust(hspace=0.7, wspace=0.5)

plt.figure(fig_fmm.number)  
plt.legend(frameon=False)
plt.xlim([0, 1])  
plt.ylim([0, 1]) 
plt.subplots_adjust(hspace=0.7, wspace=0.5)

if save_on == True:
    fig_gen.savefig("figures/im_runtime_generation_vs_resolution.pdf", format='pdf', bbox_inches='tight')   
    fig_fmm.savefig("figures/im_runtime_fmm_vs_resolution.pdf", format='pdf', bbox_inches='tight') 
    
    


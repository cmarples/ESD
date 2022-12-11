# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 17:57:40 2022

@author: Callum Marples

Plot fast marching resolution data.
Data generated in:
    - apps/run_geodesic_resolution.py
    - apps/run_geodesic_resolution_comparison.py
Note that array sizes are hardcoded in this script.

This generates Figures 4 and 5 of the Geodesic Paper.
"""

import math
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

# Convert Earth data from m to km
d_dijkstra_4[:, 1] /= 1000.0
d_fmm_4[:, 1] /= 1000.0
d_dijkstra_8[:, 1] /= 1000.0
d_fmm_8[:, 1] /= 1000.0
d_dijkstra_tri[:, 1] /= 1000.0
d_fmm_tri[:, 1] /= 1000.0
d_dijkstra_split[:, 1] /= 1000.0
d_fmm_split[:, 1] /= 1000.0
d_geo[1] /= 1000.0
d_bvm[1] /= 1000.0
d_taxi[1] /= 1000.0


m = [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]
ms = 9
lw = 1
og = [[1, 0.65, 0]]

fig_da, ax_da = plt.subplots(3, 2, figsize=set_size(fig_width, height=1.0, subplots=(3,2)))
fig_fmm, ax_fmm = plt.subplots(3, 2, figsize=set_size(fig_width, height=1.0, subplots=(3,2)))
shape_name = ["Sphere", "WGS84", "Oblate", "Prolate", "Triaxial"]
shape_label = ["(a)", "(b)", "(c)", "(d)", "(e)"]

# Plot distance data
for j in range(5): 
    
    # Plot non-wavefront distances.            
    if j == 0:
        ax_da[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [d_true[j]]*25, 'k', linestyle="dotted", linewidth=1.5)
        ax_fmm[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [d_true[j]]*25, 'k', linestyle="dotted", linewidth=1.5)
    elif j > 0:
        if j < 4:
            ax_da[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [d_geo[j]]*25, 'k', linestyle="solid", linewidth=lw)
            ax_fmm[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [d_geo[j]]*25, 'k', linestyle="solid", linewidth=lw)
        ax_da[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [d_bvm[j]]*25, 'k', linestyle="dashed", linewidth=lw)
        ax_fmm[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [d_bvm[j]]*25, 'k', linestyle="dashed", linewidth=lw)
    if j < 4:
        ax_da[m[j][0]][m[j][1]].plot(np.log10(no_vertices), [d_taxi[j]]*25, 'k', linestyle="dashdot", linewidth=lw)
        
    ax_da[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), d_dijkstra_4[:, j], c="none", edgecolor="c", marker="s", s=ms)
    ax_da[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), d_dijkstra_8[:, j], c="none", edgecolor="b", marker="o", s=ms)
    
    ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), d_fmm_4[:, j], c="c", edgecolor="c", marker="s", s=ms)
    ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), d_fmm_8[:, j], c="b", edgecolor="b", marker="o", s=ms)
    
    if j < 2: # Sphere or WGS84
        ax_da[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), d_dijkstra_tri[:, j], c="none", edgecolor=og, marker="^", s=ms)
        ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar), d_fmm_tri[:, j], c=og, edgecolor=og, marker="^", s=ms)
    else:
        # Plot icosahedral triangulation distances, avoiding repetition when splitting has no effect.
        for i in range(25):
            ax_da[m[j][0]][m[j][1]].scatter(math.log10(no_vertices_polar[i]), d_dijkstra_split[i, j], c="none", edgecolor=og, marker="^", s=ms)
            ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar[i]), d_fmm_split[i, j], c="r", edgecolor="r", marker="v", s=ms)
            if d_dijkstra_tri[i, j] != d_dijkstra_split[i, j]:
                ax_da[m[j][0]][m[j][1]].scatter(math.log10(no_vertices_polar[i]), d_dijkstra_split[i, j], c="none", edgecolor="r", marker="v", s=ms)
            if d_fmm_tri[i, j] != d_fmm_split[i, j]:    
                ax_fmm[m[j][0]][m[j][1]].scatter(np.log10(no_vertices_polar[i]), d_fmm_tri[i, j], c=og, edgecolor=og, marker="^", s=ms)
    
    
    # Set x-axis tick labels (number of vertices).  
    ax_da[m[j][0]][m[j][1]].set_xticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    ax_da[m[j][0]][m[j][1]].set_xticklabels(["$10^2$", "$10^{2.5}$", "$10^3$", "$10^{3.5}$", "$10^4$", "$10^{4.5}$", "$10^5$"])  
    ax_fmm[m[j][0]][m[j][1]].set_xticks([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    ax_fmm[m[j][0]][m[j][1]].set_xticklabels(["$10^2$", "$10^{2.5}$", "$10^3$", "$10^{3.5}$", "$10^4$", "$10^{4.5}$", "$10^5$"])    
      
    # Axis labels.
    ax_da[m[j][0]][m[j][1]].set_xlabel('Number of vertices')
    ax_fmm[m[j][0]][m[j][1]].set_xlabel('Number of vertices')
    if j == 1:
        ax_da[m[j][0]][m[j][1]].set_ylabel('Distance (km)')
        ax_fmm[m[j][0]][m[j][1]].set_ylabel('Distance (km)')
    else:
        ax_da[m[j][0]][m[j][1]].set_ylabel('Distance')
        ax_fmm[m[j][0]][m[j][1]].set_ylabel('Distance')
    
    # Label using shape names.
    ax_da[m[j][0]][m[j][1]].text(0.5,-0.42, shape_label[j], ha="center",  transform=ax_da[m[j][0]][m[j][1]].transAxes)
    ax_fmm[m[j][0]][m[j][1]].text(0.5, -0.42, shape_label[j], ha="center",  transform=ax_fmm[m[j][0]][m[j][1]].transAxes)
    ax_da[m[j][0]][m[j][1]].text(0.5, 1.05, shape_name[j], ha="center",  transform=ax_da[m[j][0]][m[j][1]].transAxes)
    ax_fmm[m[j][0]][m[j][1]].text(0.5, 1.05, shape_name[j], ha="center",  transform=ax_fmm[m[j][0]][m[j][1]].transAxes)
    
    ax_da[m[j][0]][m[j][1]].spines['top'].set_visible(False)
    ax_fmm[m[j][0]][m[j][1]].spines['top'].set_visible(False)  
    ax_da[m[j][0]][m[j][1]].spines['right'].set_visible(False)
    ax_fmm[m[j][0]][m[j][1]].spines['right'].set_visible(False)   
       
# Add legend in final ax slot.
x = -1
y = -1
j = 5

ax_da[m[j][0]][m[j][1]].scatter(x, y, c="none", edgecolor="c", marker="s", s = ms, label = "Dijkstra (polar 4)")
ax_da[m[j][0]][m[j][1]].scatter(x, y, c="none", edgecolor="b", marker="o", s = ms, label = "Dijkstra (polar 8)")
ax_da[m[j][0]][m[j][1]].scatter(x, y, c="none", edgecolor=og, marker="^", s = ms, label = "Dijkstra (icosahedral)")
ax_da[m[j][0]][m[j][1]].scatter(x, y, c="none", edgecolor="r", marker="v", s = ms, label = "Dijkstra (split icosahedral)")
ax_da[m[j][0]][m[j][1]].plot([-2, -1], [-2, -1], 'k', linestyle="dotted", label = "True distance")
ax_da[m[j][0]][m[j][1]].plot(x, y, 'k', linestyle="solid", label = "GeographicLib")
ax_da[m[j][0]][m[j][1]].plot(x, y, 'k', linestyle="dashed", label = "Boundary value method")
ax_da[m[j][0]][m[j][1]].plot(x, y, 'k', linestyle="dashdot", label = "Taxicab distance")

ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c="c", edgecolor="c", marker="s", s = ms, label = "FMM (polar 4)")
ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c="b", edgecolor="b", marker="o", s = ms, label = "FMM (polar 8)")
ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c=og, edgecolor=og, marker="^", s = ms, label = "FMM (icosahedral)")
ax_fmm[m[j][0]][m[j][1]].scatter(x, y, c="r", edgecolor="r", marker="v", s = ms, label = "FMM (split icosahedral)")
ax_fmm[m[j][0]][m[j][1]].plot([-2, -1], [-2, -1], 'k', linestyle="dotted", label = "True distance")
ax_fmm[m[j][0]][m[j][1]].plot(x, y, 'k', linestyle="solid", label = "GeographicLib")
ax_fmm[m[j][0]][m[j][1]].plot(x, y, 'k', linestyle="dashed", label = "Boundary value method")


ax_da[m[j][0]][m[j][1]].legend(loc="center left")   
ax_fmm[m[j][0]][m[j][1]].legend(loc="center left")

ax_da[m[j][0]][m[j][1]].axis('off') 
ax_fmm[m[j][0]][m[j][1]].axis('off')  

plt.figure(fig_da.number)  
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
    fig_da.savefig("figures/im_distance_dijkstra_vs_resolution.pdf", format='pdf', bbox_inches='tight')   
    fig_fmm.savefig("figures/im_distance_fmm_vs_resolution.pdf", format='pdf', bbox_inches='tight') 
    
    






      
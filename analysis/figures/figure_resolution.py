# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 19:21:01 2022

@author: Callum Marples

Plot distance as a function of resolution.

Figures 4 and 5 of the Geodesic Paper
"""

import math
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines # Used for creating a proxy artist for the legend
from figure_size import set_size

# Read sphere data
file_name_1 = "../data/distance_vs_resolution_sphere.csv"
L = [0] * 10
L_no = 0
with open(file_name_1) as f:
    
    csv_reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in csv_reader:
        if L_no < 10:
            L[L_no] = row
            L_no += 1

d_1s, d_2s, d_3s, d_4s = L[0], L[1], L[2], L[3]
t_1s, t_2s, t_3s, t_4s = L[4], L[5], L[6], L[7]
a_s, b_s, c_s = L[8][0], L[8][1], L[8][2]
d_true_s = L[9]

# Read prolate data
file_name_2 = "../data/distance_vs_resolution_prolate.csv"
L = [0] * 10
L_no = 0
with open(file_name_2) as f:
    
    csv_reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in csv_reader:
        if L_no < 10:
            L[L_no] = row
            L_no += 1

d_1p, d_2p, d_3p, d_4p = L[0], L[1], L[2], L[3]
t_1p, t_2p, t_3p, t_4p = L[4], L[5], L[6], L[7]
a_p, b_p, c_p = L[8][0], L[8][1], L[8][2]
d_geo_lib_p = L[9]

# Read triaxial data
file_name_3 = "../data/distance_vs_resolution_triaxial.csv"
L = [0] * 10
L_no = 0
with open(file_name_3) as f:
    
    csv_reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in csv_reader:
        if L_no < 10:
            L[L_no] = row
            L_no += 1

d_1t, d_2t, d_3t, d_4t = L[0], L[1], L[2], L[3]
t_1t, t_2t, t_3t, t_4t = L[4], L[5], L[6], L[7]
a_t, b_t, c_t = L[8][0], L[8][1], L[8][2]
d_bvm_t = L[9]

# Read Earth data
file_name_4 = "../data/distance_vs_resolution_earth.csv"
L = [0] * 10
L_no = 0
with open(file_name_4) as f:
    
    csv_reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in csv_reader:
        if L_no < 10:
            L[L_no] = row
            L_no += 1

d_1e, d_2e, d_3e, d_4e = L[0], L[1], L[2], L[3]
t_1e, t_2e, t_3e, t_4e = L[4], L[5], L[6], L[7]
a_e, b_e, c_e = L[8][0], L[8][1], L[8][2]
d_geo_lib_e = L[9]

# Dependent variable (number of phi values, with n_theta = n_phi+1)
n_div = np.arange(25, 525, step=25)
m = len(n_div)
s1, s2 = 26, 20 # Number of values in 'x' and 'y'arrays
n_div_ticks = np.arange(100, 600, step=100)

# Set figure parameters
plt.style.use('tex')
width = 345.0
fs = set_size(width, 1.0, 1.0, [3, 2])
# Initialise figure objects
fig = plt.figure(figsize=(fs[0], fs[1])) 
gs = plt.GridSpec(3, 2)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])
ax4 = fig.add_subplot(gs[3])
ax5 = fig.add_subplot(gs[4])
ax6 = fig.add_subplot(gs[5])
fig.subplots_adjust(wspace=0.5, hspace=0.25)

# Plot the distance figures
# Sphere
line1 = ax1.scatter(n_div, d_1s, c="none", edgecolor="c", marker="o", s = s1)
line2 = ax1.scatter(n_div, d_2s, c="none", edgecolor="m", marker="s", s = s2)
line3 = ax1.scatter(n_div, d_3s, c="c", edgecolor="c", marker="o", s = s1) 
line4 = ax1.scatter(n_div, d_4s, c="m", edgecolor="m", marker="s", s = s2)
ax1.plot([25, 500], [d_true_s, d_true_s], 'k')
ax1.set_xticks(n_div_ticks) 
ax1.set_xlabel('$n_{\phi}$')
ax1.set_ylabel('Distance')
ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)

# Earth
ax2.scatter(n_div, d_1e, c="none", edgecolor="c", marker="o", s = s1)
ax2.scatter(n_div, d_2e, c="none", edgecolor="m", marker="s", s = s2)
ax2.scatter(n_div, d_3e, c="c", edgecolor="c", marker="o", s = s1) 
ax2.scatter(n_div, d_4e, c="m", edgecolor="m", marker="s", s = s2)
line6 = ax2.plot([25, 500], [d_geo_lib_e, d_geo_lib_e], 'k--')
ax2.set_xticks(n_div_ticks) 
ax2.set_xlabel('$n_{\phi}$')
ax2.set_ylabel('Distance (km)')
ax2.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)

# Prolate
ax3.scatter(n_div, d_1p, c="none", edgecolor="c", marker="o", s = s1)
ax3.scatter(n_div, d_2p, c="none", edgecolor="m", marker="s", s = s2)
ax3.scatter(n_div, d_3p, c="c", edgecolor="c", marker="o", s = s1) 
ax3.scatter(n_div, d_4p, c="m", edgecolor="m", marker="s", s = s2)
ax3.plot([25, 500], [d_geo_lib_p, d_geo_lib_p], 'k--')
ax3.set_xticks(n_div_ticks) 
ax3.set_xlabel('$n_{\phi}$')
ax3.set_ylabel('Distance')
ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)

# Oblate
ax4.scatter(n_div, d_1p, c="none", edgecolor="c", marker="o", label='First Order', s = s1)
ax4.scatter(n_div, d_2p, c="none", edgecolor="m", marker="s", label='Second Order', s = s2)
ax4.scatter(n_div, d_3p, c="c", edgecolor="c", marker="o", label='First Order with Refinement', s = s1) 
ax4.scatter(n_div, d_4p, c="m", edgecolor="m", marker="s", label='Second Order with Refinement', s = s2)
ax4.plot([25, 500], [d_geo_lib_p, d_geo_lib_p], 'k--', label='GeographicLib')
ax4.set_xticks(n_div_ticks) 
ax4.set_xlabel('$n_{\phi}$')
ax4.set_ylabel('Distance')
ax4.spines["right"].set_visible(False)
ax4.spines["top"].set_visible(False)

# Triaxial
ax5.scatter(n_div, d_1t, c="none", edgecolor="c", marker="o", label='First Order', s = s1)
ax5.scatter(n_div, d_2t, c="none", edgecolor="m", marker="s", label='Second Order', s = s2)
ax5.scatter(n_div, d_3t, c="c", edgecolor="c", marker="o", label='First Order with Refinement', s = s1) 
ax5.scatter(n_div, d_4t, c="m", edgecolor="m", marker="s", label='Second Order with Refinement', s = s2)
line7 = ax5.plot([25, 500], [d_bvm_t, d_bvm_t], 'k:', label='Boundary value method')
ax5.set_xticks(n_div_ticks) 
ax5.set_xlabel('$n_{\phi}$')
ax5.set_ylabel('Distance')
ax5.spines["right"].set_visible(False)
ax5.spines["top"].set_visible(False)

# Legend
line5 = mlines.Line2D([], [], color='black')
line6 = mlines.Line2D([], [], color='black', linestyle='--')
line7 = mlines.Line2D([], [], color='black', linestyle=':')
ax6.legend([line1, line2, line3, line4, line5, line6, line7], ['First order',
                                                               'Second Order',
                                                               'First Order with Refinement',
                                                               'Second Order with Refinement',
                                                               'True distance',
                                                               'GeographicLib',
                                                               'Boundary value method'],
           loc="center")
ax6.xaxis.set_visible(False)
ax6.yaxis.set_visible(False)
for v in ax6.spines.values():
    v.set_visible(False)
ax1.text(0.7, 0.9, 'Unit sphere', horizontalalignment='center',
     verticalalignment='center', transform=ax1.transAxes)
ax2.text(0.7, 0.9, 'Earth (WGS84)', horizontalalignment='center',
     verticalalignment='center', transform=ax2.transAxes)
ax3.text(0.7, 0.8, 'Prolate (1,1,2)', horizontalalignment='center',
     verticalalignment='center', transform=ax3.transAxes)
ax4.text(0.7, 0.8, 'Oblate (2,2,1)', horizontalalignment='center',
     verticalalignment='center', transform=ax4.transAxes)
ax5.text(0.7, 0.9, 'Triaxial (3,2,1)', horizontalalignment='center',
     verticalalignment='center', transform=ax5.transAxes)

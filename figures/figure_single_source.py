# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 21:56:48 2022

@author: Callum Marples

Plot single source fast marching data.
Data generated in:
    - apps/run_geodesic_resolution.py
    - apps/run_geodesic_resolution_comparison.py
Note that array sizes are hardcoded in this script.

These are Figures 8 and 9 of the Geodesic Paper.
"""

import math
import numpy as np
import csv
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import ListedColormap
from matplotlib import colors
from figure_size import set_size

print("Plot single source distance difference maps")

os.chdir("..")

save_on = True
plt.style.use('tex')
fig_width = 345.0

fig1 = plt.figure(figsize=set_size(fig_width, height=1.0))
plt.subplots_adjust(wspace=0.3, hspace=0.3)
fig2 = plt.figure(figsize=set_size(fig_width, height=1.0))


a1, b1, c1 = 1.0, 1.0, 1.0 # Sphere axes
a2, b2, c2 = 3.0, 2.0, 1.0 # Ellipsoid axes
r = (a2*b2*c2) ** (1/3)
a2 /= r
b2 /= r
c2 /= r

start_point = [90.0, 0.0]

### Read data
n = 308 # Number of endpoints
end_point = [0] * n
d_sphere_polar = [0] * n
d_sphere_ico = [0] * n
d_sphere_true = [0] * n
d_triaxial_polar = [0] * n
d_triaxial_ico = [0] * n
d_triaxial_bvm = [0] * n
t_true = [0] * n
t_bvm = [0] * n

# Read data from sphere file
file_name = "data/geodesics/single_source_sphere_true.txt"
r_no = 0
with open(file_name) as f:
    csv_reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in csv_reader:
        if r_no < n:
            end_point[r_no] = [row[0], row[1]]
            d_sphere_true[r_no] = row[2]
            t_true[r_no] = row[3]
            r_no += 1

# Read data from triaxial file
file_name = "data/geodesics/single_source_ellipsoid_bvm.txt"
r_no = 0
with open(file_name) as f:
    csv_reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in csv_reader:
        if r_no < n:
            d_triaxial_bvm[r_no] = row[2]
            t_bvm[r_no] = row[3]
            r_no += 1
            
# Read data from FMM file
file_name = "data/geodesics/single_source_fast_marching.txt"
r_no = 0
with open(file_name) as f:
    csv_reader = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in csv_reader:
        if r_no < n:
            d_sphere_polar[r_no] = row[0]
            d_sphere_ico[r_no] = row[1]
            d_triaxial_polar[r_no] = row[2]
            d_triaxial_ico[r_no] = row[3]
        elif r_no == n:
            t_sphere_polar_gen = row[0]
            t_sphere_ico_gen = row[1]
            t_triaxial_polar_gen = row[2]
            t_triaxial_ico_gen = row[3]
        elif r_no == n+1:
            t_sphere_polar_fmm = row[0]
            t_sphere_ico_fmm = row[1]
            t_triaxial_polar_fmm = row[2]
            t_triaxial_ico_fmm = row[3]
        r_no += 1 
            
### Deal with gaps in the triaxial boundary value method data 
D = dict()
for i in range(n):
    D[(end_point[i][0], end_point[i][1])] = d_triaxial_bvm[i]
                
m = 0
d_gaps = []
for i in range(n):
    if d_triaxial_bvm[i] == -1.0:
        d_gaps.append(i)
        #print(end_angles[i])
        th = end_point[i][0]
        ph = end_point[i][1]
        if ph == 180 or ph == 0:
            d_triaxial_bvm[i] = D[(180 - end_point[i][0], end_point[i][1])]
            D[(end_point[i][0], end_point[i][1])] = d_triaxial_bvm[i]
        elif ph != 0:
            d_triaxial_bvm[i] = D[(end_point[i][0], 360 - end_point[i][1])]
            D[(end_point[i][0], end_point[i][1])] = d_triaxial_bvm[i]
        m += 1    
    

### Difference map
def difference_map(d_true, d_fmm):
    
    diff = [0] * n
    for i in range(n):
        diff[i] = d_fmm[i] - d_true[i]
        if d_true[i] == -1.0:
            diff[i] = math.nan
        
    d_mat = np.zeros((19, 18))
    d_mat[1:-1] = np.reshape(diff[1:-1], (17, 18))
        
    d2 = np.zeros((19, 19))
    d2[:,:-1] = d_mat
    d2[:,-1] = d2[:,0]
    
    d2[0] = diff[0]
    d2[-1] = diff[-1] 
    
    return d2, d_mat

### Colour image plotting routines
def plot_map_2D(ax, d2, cmap, norm):
    ax.imshow(d2, cmap=cmap, norm=norm) 
    ax.plot(0, 8.9, color=[0,0,0], marker='*', markersize=5)
    ax.plot(18, 8.9, color=[0,0,0], marker='*', markersize=5)
    xx = np.arange(0, 19, 3)
    ax.set_xticks(xx)
    ax.set_xticklabels(["0", "60", "120", "180", "240", "300", "360"])
    yy = np.arange(0, 19, 3)
    ax.set_yticks(yy)
    ax.set_yticklabels(["0", "30", "60", "90", "120", "150", "180"])
    ax.set_xlabel('$\phi$ (degrees)')
    ax.set_ylabel('$\\theta$ (degrees)')

def plot_map_2D_triaxial(ax, d2, cmap, norm, d_gaps, marker, s):
    ax.imshow(d2, cmap=cmap, norm=norm) 
    ax.plot(0, 8.9, color=[0,0,0], marker='*', markersize=5)
    ax.plot(18, 8.9, color=[0,0,0], marker='*', markersize=5)
    xx = np.arange(0, 19, 3)
    ax.set_xticks(xx)
    ax.set_xticklabels(["0", "60", "120", "180", "240", "300", "360"])
    yy = np.arange(0, 19, 3)
    ax.set_yticks(yy)
    ax.set_yticklabels(["0", "30", "60", "90", "120", "150", "180"])
    for i in d_gaps:
        yg, xg = end_point[i]
        yg = yg/10
        xg = xg/20
        ax.plot(xg, yg, color='k', marker=marker, markersize=s)
        if yg == 0 or yg == 18:
            for j in range(19):
                ax.plot(j, yg, color='k', marker=marker, markersize=s)
        elif xg == 0:
            ax.plot(18, yg, color='k', marker=marker, markersize=s)
                
        
def plot_map_ellipsoid(ax, d_mat, d2, a, b, c, cmap, norm):
    u, v = np.mgrid[0:np.pi:19j, 0:2*np.pi:19j]
    v -= 10*np.pi/180
    for (i, j), value in np.ndenumerate(d_mat):
        if d2[i][j] == -1.0:
            d_mat[i][j] = d_mat[9][9]
    strength = d_mat
    x = a * np.sin(u) * np.cos(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(u)
    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cmap, linewidth=0, antialiased=True,
                           facecolors=cmap(norm(strength)))
    ax.scatter(a+0.1, 0, 0, color='k', marker='*', s=12)
    ax.set_box_aspect((1,1,1))
    ax.set_axis_off()


### Sphere figure
d2_sphere_polar, d_mat_sphere_polar = difference_map(d_sphere_true, d_sphere_polar)
d2_sphere_ico, d_mat_sphere_ico = difference_map(d_sphere_true, d_sphere_ico)

cmap = mpl.cm.RdYlBu

ax1c = fig1.add_axes([0.56, 0.545, 0.015, 0.335])
ax2c = fig1.add_axes([0.56, 0.11, 0.015, 0.335])
lim = 0.012

norm = mpl.colors.Normalize(vmin=-lim, vmax=lim)
cb1 = mpl.colorbar.ColorbarBase(ax1c, cmap=cmap, norm=norm)
cb2 = mpl.colorbar.ColorbarBase(ax2c, cmap=cmap, norm=norm)
         
ax1 = fig1.add_subplot(2, 3, (1,2))
plot_map_2D(ax1, d2_sphere_polar, cmap, norm)
ax1.set_xlabel('$\phi$ (degrees)')
ax1.set_ylabel('$\\theta$ (degrees)')

ax4 = fig1.add_subplot(2, 3, (4,5))
plot_map_2D(ax4, d2_sphere_ico, cmap, norm)
ax4.set_xlabel('$\phi$ (degrees)')
ax4.set_ylabel('$\\theta$ (degrees)')

ax2 = fig1.add_subplot(4, 3, 3, projection='3d')
plot_map_ellipsoid(ax2, d_mat_sphere_polar, d2_sphere_polar, a1, b1, c1, cmap, norm)
ax2.view_init(elev=15, azim=20)

ax3 = fig1.add_subplot(4, 3, 6, projection='3d')
plot_map_ellipsoid(ax3, d_mat_sphere_polar, d2_sphere_polar, a1, b1, c1, cmap, norm)
ax3.view_init(elev=15, azim=200)

ax5 = fig1.add_subplot(4, 3, 9, projection='3d')
plot_map_ellipsoid(ax5, d_mat_sphere_ico, d2_sphere_ico, a1, b1, c1, cmap, norm)
ax5.view_init(elev=15, azim=20)
ax6 = fig1.add_subplot(4, 3, 12, projection='3d')
plot_map_ellipsoid(ax6, d_mat_sphere_ico, d2_sphere_ico, a1, b1, c1, cmap, norm)
ax6.view_init(elev=15, azim=200)

ax1.text(-0.5, 0.5, "(a)", ha="center",  transform=ax1.transAxes)
ax4.text(-0.5, 0.5, "(b)", ha="center",  transform=ax4.transAxes)

### Ellipsoid figure
d2_triaxial_polar, d_mat_triaxial_polar = difference_map(d_triaxial_bvm, d_triaxial_polar)
d2_triaxial_ico, d_mat_triaxial_ico = difference_map(d_triaxial_bvm, d_triaxial_ico)

upper_lim = math.pi*b2*c2/a2
d_bvm_valid = [0] * 308
for i in range(308):
    if d_triaxial_bvm[i] > upper_lim or d_triaxial_bvm[i] == -1.0:
        d_bvm_valid[i] = math.nan
    else:
        d_bvm_valid[i] = d_triaxial_bvm[i]

d_gaps2 = []
for i in range(n):
    if math.isnan(d_bvm_valid[i]):
        d_gaps2.append(i)
        
d2_triaxial_polar2, dmat1 = difference_map(d_bvm_valid, d_triaxial_polar)
d2_triaxial_ico2, dmat2 = difference_map(d_bvm_valid, d_triaxial_ico)


bx1c = fig2.add_axes([0.42, 0.567, 0.015, 0.276])
bx2c = fig2.add_axes([0.92, 0.567, 0.015, 0.276])
bx3c = fig2.add_axes([0.42, 0.147, 0.015, 0.276])
bx4c = fig2.add_axes([0.92, 0.147, 0.015, 0.276])

lim1 = 0.7
lim2 = 0.04
norm1 = mpl.colors.Normalize(vmin=-lim1, vmax=lim1)
norm2 = mpl.colors.Normalize(vmin=-lim2, vmax=lim2)
cb1 = mpl.colorbar.ColorbarBase(bx1c, cmap=cmap, norm=norm1)
cb2 = mpl.colorbar.ColorbarBase(bx2c, cmap=cmap, norm=norm2)
cb3 = mpl.colorbar.ColorbarBase(bx3c, cmap=cmap, norm=norm1)
cb4 = mpl.colorbar.ColorbarBase(bx4c, cmap=cmap, norm=norm2)

bx1 = fig2.add_subplot(2, 2, 1)
plot_map_2D_triaxial(bx1, d2_triaxial_polar, cmap, norm1, d_gaps, 'x', 2.2)
bx1.set_ylabel('$\\theta$ (degrees)')
bx1.set_title('(a)', fontsize=10)

bx2 = fig2.add_subplot(2, 2, 2)
plot_map_2D_triaxial(bx2, d2_triaxial_polar2, cmap, norm2, d_gaps2, '_', 2)
bx2.set_title('(b)', fontsize=10)

bx3 = fig2.add_subplot(2, 2, 3)
plot_map_2D_triaxial(bx3, d2_triaxial_ico, cmap, norm1, d_gaps, 'x', 2.2)
bx3.set_xlabel('$\phi$ (degrees)')
bx3.set_ylabel('$\\theta$ (degrees)')
bx3.set_title('(c)', fontsize=10)

bx4 = fig2.add_subplot(2, 2, 4)
plot_map_2D_triaxial(bx4, d2_triaxial_ico2, cmap, norm2, d_gaps2, '_', 2)
bx4.plot(0.99, 15, color='w', marker='x', markersize=3)
bx4.plot(17.1, 15, color='w', marker='x', markersize=3)
bx4.set_xlabel('$\phi$ (degrees)')
bx4.set_title('(d)', fontsize=10)

plt.subplots_adjust(wspace=0.8, hspace=0.2)

if save_on:
    fig1.savefig('figures/im_single_source_sphere.pdf', format='pdf', bbox_inches='tight')
    fig2.savefig('figures/im_single_source_triaxial.pdf', format='pdf', bbox_inches='tight')
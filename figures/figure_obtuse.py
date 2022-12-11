# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 23:48:04 2022

@author: Callum Marples

Draw the mesh of a triaxial ellipsoid using each of the two methods; polar and
icosahedral. Colour each triangle/quadrilateral according to its largest angle.

This is Figure 3 of the Geodesic Paper.
"""

import leod
import numpy as np
import math
import matplotlib.pyplot as plt
from figure_size import set_size
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import ListedColormap

### Initialise figure
save_on = True
plt.style.use('tex')
fig_width = 345.0

fig, axes = plt.subplots(1, 3, figsize=set_size(fig_width, height=1.0, subplots=(1, 3)))
ax = axes.flat

### Define shape
a = 3.0
b = 2.0
c = 1.0
shape = leod.ellipsoid_shape.EllipsoidShape(a, b, c)

rad2deg = 180.0 / math.pi # Conversion factor

### Polar mesh
no_theta = 9
no_phi = 16
vertex, polar_grid = leod.fmm_polar_graph.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=False)
no_vertices = len(vertex)
leod.fmm_precalculation.precalculate_grid(vertex)
[no_obtuse, max_angle] = leod.fmm_precalculation.find_obtuse_angles(vertex)

ax[0].remove()
ax[0] = fig.add_subplot(1, 3, 1, projection='3d')

for i in range(no_vertices):
    
    if i == 0 or i == len(vertex)-1:
        for j in vertex[i].neighbour.keys():

            for k in vertex[i].neighbour[j].face:

                x = [vertex[i].carts[0], vertex[j].carts[0], vertex[k].carts[0]]
                y = [vertex[i].carts[1], vertex[j].carts[1], vertex[k].carts[1]]
                z = [vertex[i].carts[2], vertex[j].carts[2], vertex[k].carts[2]]
                verts = [list(zip(x, y, z))]
                
                # Find largest angle in face
                cos_alpha_i = vertex[i].neighbour[j].face_angle[k]
                cos_alpha_j = vertex[j].neighbour[i].face_angle[k]
                cos_alpha_k = vertex[k].neighbour[i].face_angle[j]
                cos_alpha = cos_alpha_i
                if cos_alpha_j < cos_alpha:
                    cos_alpha = cos_alpha_j
                if cos_alpha_k < cos_alpha:
                    cos_alpha = cos_alpha_k
                alpha_max = math.acos(cos_alpha) * rad2deg
                    
                if alpha_max <= 90.0:
                    face_col = [1, 1, (alpha_max-60.0)/30.0]
                    srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                    ax[0].add_collection3d(srf)
                    
                else:
                    face_col = [1 - (alpha_max-75.0)/90.0, 1 - (alpha_max-75.0)/90.0, 1]
                    srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                    ax[0].add_collection3d(srf)
    
    else:
        th = leod.fmm_polar_graph.get_theta_index(i, no_vertices, no_theta, no_phi)
        if th < no_theta - 2:
            ph = leod.fmm_polar_graph.get_phi_index(i, th, no_phi)
            ph_1 = ph + 1
            if ph_1 == no_phi:
                ph_1 -= no_phi
            th_1 = th + 1
            j = leod.fmm_polar_graph.get_vertex_index(th, ph_1, no_vertices, no_theta, no_phi)
            k = leod.fmm_polar_graph.get_vertex_index(th_1, ph, no_vertices, no_theta, no_phi)
            l = leod.fmm_polar_graph.get_vertex_index(th_1, ph_1, no_vertices, no_theta, no_phi)
            
            cos_alpha_i = vertex[i].neighbour[j].face_angle[k]
            cos_alpha_j = vertex[j].neighbour[i].face_angle[l]
            cos_alpha_k = vertex[k].neighbour[i].face_angle[l]
            cos_alpha_l = vertex[l].neighbour[j].face_angle[k]
            cos_alpha = cos_alpha_i
            if cos_alpha_j < cos_alpha:
                cos_alpha = cos_alpha_j
            if cos_alpha_k < cos_alpha:
                cos_alpha = cos_alpha_k
            if cos_alpha_l < cos_alpha:
                cos_alpha = cos_alpha_l
            alpha_max = math.acos(cos_alpha) * rad2deg
                    
            x = [vertex[i].carts[0], vertex[j].carts[0], vertex[l].carts[0], vertex[k].carts[0]]
            y = [vertex[i].carts[1], vertex[j].carts[1], vertex[l].carts[1], vertex[k].carts[1]]
            z = [vertex[i].carts[2], vertex[j].carts[2], vertex[l].carts[2], vertex[k].carts[2]]
            verts = [list(zip(x, y, z))]
            
            if alpha_max <= 90.0:
                face_col = [1, 1, (alpha_max-60.0)/30.0]
                srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                ax[0].add_collection3d(srf)
                
            else:
                face_col = [1 - (alpha_max-75.0)/90.0, 1 - (alpha_max-75.0)/90.0, 1]
                #face_col = [1 - (alpha_max-90.0)/45.0, 1 - (alpha_max-90.0)/45.0, 1]
                srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                ax[0].add_collection3d(srf)

limit = 1.6
ax[0].set_xlim([-limit, limit])
ax[0].set_ylim([-limit, limit])
ax[0].set_zlim([-limit, limit])
ax[0].set_title('(a)', y=-0.01, fontsize=10)
ax[0].set_axis_off()

### Colour bar

ax[1].remove()
ax[1] = fig.add_axes([0.47, 0.2, 0.02, 0.6])

cmap = mpl.cm.cool


n = 75
x = np.linspace(0.0, 1.0, n)
cm = np.zeros([n, 4])
for i in range(n):
    if i <= 29:
        cm[i][0] = 1.0
        cm[i][1] = 1.0
        cm[i][2] = (i+60 - 60.0) / 30.0
    else:
        cm[i][0] = 1.0 - (i+60-75.0)/90.0
        cm[i][1] = 1.0 - (i+60-75.0)/90.0
        cm[i][2] = 1.0
    cm[i][3] = 1.0

newcmp = ListedColormap(cm)
norm = mpl.colors.Normalize(vmin=60, vmax=135)

cb1 = mpl.colorbar.ColorbarBase(ax[1], cmap=newcmp, norm=norm, orientation='vertical',
                                ticks=[60, 75, 90, 105, 120, 135] )
cb1.set_label(r'$\psi_{\mathrm{max}}$ (degrees)')






### Icosahedral triangulation

# Generate triangulation for the unit sphere.
n = 3  # number of triangular divisions.
vertex2 = leod.triangulation_sphere.triangulate_sphere(1.0, n)

# Scale vertices to the surface of the ellipsoid.
for i in range(len(vertex2)):
    vertex2[i].carts[0] *= shape.a_axis
    vertex2[i].carts[1] *= shape.b_axis
    vertex2[i].carts[2] *= shape.c_axis

# Calculate face angles.
leod.fmm_precalculation.precalculate_grid(vertex2)
[no_obtuse2, max_angle2] = leod.fmm_precalculation.find_obtuse_angles(vertex2)


# Plot faces of the triangulation
ax[2].remove()
ax[2] = fig.add_subplot(1, 3, 3, projection='3d')



for i in range(len(vertex2)):
    for j in vertex2[i].neighbour.keys():
        if j > i:
            for k in vertex2[i].neighbour[j].face:
                if k > j:
                    # Find largest angle in face
                    cos_alpha_i = vertex2[i].neighbour[j].face_angle[k]
                    cos_alpha_j = vertex2[j].neighbour[i].face_angle[k]
                    cos_alpha_k = vertex2[k].neighbour[i].face_angle[j]
                    cos_alpha = cos_alpha_i
                    if cos_alpha_j < cos_alpha:
                        cos_alpha = cos_alpha_j
                    if cos_alpha_k < cos_alpha:
                        cos_alpha = cos_alpha_k
                    alpha_max = math.acos(cos_alpha) * rad2deg
                    
                    # Get face vertices in the appropriate format
                    x = [vertex2[i].carts[0], vertex2[j].carts[0], vertex2[k].carts[0]]
                    y = [vertex2[i].carts[1], vertex2[j].carts[1], vertex2[k].carts[1]]
                    z = [vertex2[i].carts[2], vertex2[j].carts[2], vertex2[k].carts[2]]
                    verts = [list(zip(x, y, z))]
                    
                    # Determine face colour using the angle and plot face
                    if alpha_max <= 90.0:
                        face_col = [1, 1, (alpha_max-60.0)/30.0]
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                        ax[2].add_collection3d(srf)
                        
                    else:
                        #face_col = [1, 1 - (alpha_max-90.0)/90.0, 1]
                        face_col = [1 - (alpha_max-75.0)/90.0, 1 - (alpha_max-75.0)/90.0, 1]
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                        ax[2].add_collection3d(srf)
                    
limit = 1.6
ax[2].set_xlim([-limit, limit])
ax[2].set_ylim([-limit, limit])
ax[2].set_zlim([-limit, limit])
ax[2].set_title('(b)', y=-0.01, fontsize=10)
ax[2].set_axis_off()

fig.subplots_adjust(hspace=0.1, wspace=0.1)

if save_on:
    fig.savefig('im_obtuse_angles.pdf', format='pdf', bbox_inches='tight') 
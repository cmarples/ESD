# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:15:35 2022

@author: Callum Marples

Generate figure showing the definition of the theta-phi mesh.

This is Figure 1 of the Geodesic Paper.
"""

import leod
import numpy as np
import math
import matplotlib.pyplot as plt
from figure_size import set_size
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

save_on = True
plt.style.use('tex')
fig_width = 345.0

n_th, n_ph = 5, 4

N = (n_th-2)*n_ph + 2

index_dict = {1: (0, 0),
              2: (1, 0),
              3: (1, 1),
              4: (1, 2),
              5: (1, 3),
              6: (2, 0),
              7: (2, 1),
              8: (2, 2),
              9: (2, 3),
              10: (3, 0),
              11: (3, 1),
              12: (3, 2),
              13: (3, 3),
              14: (4, 0)}


th = np.linspace(0, 180, n_th)
ph0 = np.linspace(0, 360, n_ph+1)

ph = ph0 + 180/n_ph

fig, axes = plt.subplots(1, 2, figsize=set_size(fig_width, height=1.0, subplots=(1, 2)))

ax = axes.flat

#fig = plt.figure(figsize=set_size(fig_width, height=1.0))
#fig = plt.figure(figsize=plt.figaspect(0.5))
#ax = fig.add_subplot(1, 2, 1)

# Pole connections
for j in range(n_ph):
    ax[0].plot([180, ph[j]], [0, th[1]], 'k')
    ax[0].plot([180, ph[j]], [180, th[n_th-2]], 'k')
    
k = 1
for i in range(1, n_th-1, 1):
    for j in range(n_ph):
        k += 1
        if j == n_ph-1:
            ax[0].plot([ph[j], ph[j]+0.5*ph[1]], [th[i], th[i]], 'k', linestyle='dashed')
            if i < n_th/2:
                ax[0].plot([ph[j], ph[j]+0.5*ph[1]], [th[i], th[i]+th[1]], 'k', linestyle='dashdot')
                ax[0].plot([ph[j], ph[j]+0.5*ph[1]], [th[i]+th[1], th[i]], 'k', linestyle='dashdot')
        else:
            if j == 0:
                ax[0].plot([ph[j], ph[j]-0.5*ph[1]], [th[i], th[i]], 'k', linestyle='dashed')
                if i < n_th/2:
                    ax[0].plot([ph[j], ph[j]-0.5*ph[1]], [th[i], th[i]+th[1]], 'k', linestyle='dashdot')
                    ax[0].plot([ph[j], ph[j]-0.5*ph[1]], [th[i]+th[1], th[i]], 'k', linestyle='dashdot')
            ax[0].plot([ph[j], ph[j+1]], [th[i], th[i]], 'k', linestyle='solid')
            
            if i > 0 and i < n_th-2:
                ax[0].plot([ph[j], ph[j+1]], [th[i], th[i+1]], 'k', linestyle='dotted')
                ax[0].plot([ph[j], ph[j+1]], [th[i+1], th[i]], 'k', linestyle='dotted')
            
        if i < n_th-2:
                ax[0].plot([ph[j], ph[j]], [th[i], th[i+1]], 'k', linestyle='solid')
                
        ax[0].plot(ph[j], th[i], 'b.', markersize=10)

ax[0].plot(180, 0, 'b.', markersize=10)
ax[0].plot(180, 180, 'b.', markersize=10)
ax[0].text(80, 200,  r'South pole ($\theta = \pi$)')  
ax[0].text(80, -10, r'North pole ($\theta = 0$)') 
  

ax[0].text(10, 25, r"$\phi = 0$")

ax[0].invert_yaxis()

# Remove axes
ax[0].axis('off')
ax[0].set_title('(a)', y=-0.2, fontsize=10)

#plt.axis('off')
#plt.show()

# Generate sphere plot

#fig = plt.figure(figsize=(fs[0], fs[0]))
ax[1].remove()
ax[1] = fig.add_subplot(1, 2, 2, projection='3d')

a, b, c = 3.0, 2.0, 1.0
thS = np.linspace(0, math.pi, 100)
phS = np.linspace(0, 2*math.pi, 100)
TH, PH = np.meshgrid(thS, phS)

X = a * np.sin(TH) * np.cos(PH)
Y = b * np.sin(TH) * np.sin(PH)
Z = c * np.cos(TH)

#ax[1].plot_surface(X, Y, Z, color='k', linewidth=0, alpha=0.1)

shape = leod.ellipsoid_shape.EllipsoidShape(a, b, c)
no_theta = 15
no_phi = 16
vertex, polar_grid = leod.fmm_polar_graph.generate_polar_graph(shape, no_theta, no_phi, is_connect_8=True)
no_vertices = len(vertex)

for i in range(no_vertices):
    
    if i == 0 or i == len(vertex)-1:
        for j in vertex[i].neighbour.keys():

            for k in vertex[i].neighbour[j].face:

                x = [vertex[i].carts[0], vertex[j].carts[0], vertex[k].carts[0]]
                y = [vertex[i].carts[1], vertex[j].carts[1], vertex[k].carts[1]]
                z = [vertex[i].carts[2], vertex[j].carts[2], vertex[k].carts[2]]
                verts = [list(zip(x, y, z))]
                
                srf = Poly3DCollection(verts, alpha=1.0, facecolor=[0.9, 0.9, 0.9], edgecolor='k')
                ax[1].add_collection3d(srf)
    
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
            
            x = [vertex[i].carts[0], vertex[j].carts[0], vertex[l].carts[0], vertex[k].carts[0]]
            y = [vertex[i].carts[1], vertex[j].carts[1], vertex[l].carts[1], vertex[k].carts[1]]
            z = [vertex[i].carts[2], vertex[j].carts[2], vertex[l].carts[2], vertex[k].carts[2]]
            verts = [list(zip(x, y, z))]
            
            srf = Poly3DCollection(verts, alpha=1.0, facecolor=[0.9, 0.9, 0.9], edgecolor='k')
            ax[1].add_collection3d(srf)
            
            
    
    ax[1].scatter(vertex[i].carts[0], vertex[i].carts[1], vertex[i].carts[2], color='b', s=6)
    
        
        
ax[1].view_init(elev=40., azim=61)
# Set axis limits - make all of them the same to get same scaling
#m = max([a, b, c])
m = 1.6
ax[1].set_xlim3d(-m, m)
ax[1].set_ylim3d(-m, m)
ax[1].set_zlim3d(-m, m)

# Remove axes
ax[1].axis('off')
ax[1].set_title('(b)', y=-0.2, fontsize=10)

fig.subplots_adjust(hspace=0.5,wspace=0.1)

if save_on:
    fig.savefig('im_polar_mesh.pdf', format='pdf', bbox_inches='tight')
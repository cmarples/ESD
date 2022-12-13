# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 12:32:37 2022

@author: Callum Marples

Generate figure showing the definition of the icosahedral mesh.

This is Figure 2 of the Geodesic Paper.
"""

import leod
import numpy as np
import math
import matplotlib.pyplot as plt
from figure_size import set_size
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

### Initialise figure
save_on = True
plt.style.use('tex')
fig_width = 345.0

fig, axes = plt.subplots(2, 2, figsize=set_size(fig_width, height=1.0, subplots=(2, 2)))
ax = axes.flat

ylw = [1.0, 1.0, 0.6]
blu = [0.8, 0.8, 1.0]

### Plot icosahedron

# Icosahedron coordinates
t = 0.5 * (1.0 + math.sqrt(5.0)) # Golden ratio
den = math.sqrt(t*t + 1.0)     # Divide by this to get unit circumradius

ico_vertex = [ [0, 1, t], [0, -1, t], [0, 1, -t], [0, -1, -t],
               [1, t, 0], [-1, t, 0], [1, -t, 0], [-1, -t, 0],
               [t, 0, 1], [-t, 0, 1], [t, 0, -1], [-t, 0, -1] ]
ico_vertex = np.array(ico_vertex)

# Determine neighbouring vertices (edge length = 2)
neigh = []
for i in range(12):
    neigh.append([-1] * 5)
    neighbour_no = 0
    for j in range(12):
        if i != j:
            d = np.linalg.norm(ico_vertex[i]-ico_vertex[j])
            if math.fabs(d - 2.0) < 1e-9:
                neigh[i][neighbour_no] = j
                neighbour_no += 1

# Determine faces by considering each trio in turn
ico_face = []
for i in range(10):
    for j in range(i, 11):
        for k in range(j, 12):
            if j in neigh[i] and k in neigh[i] and i in neigh[j] and k in neigh[j] and i in neigh[k] and j in neigh[k]:
                ico_face.append([i, j, k])

ylw_index = 6 # Make this face yellow
triangles = []
for i in range(len(ico_face)):
    if i != ylw_index:
        triangles.append(( (ico_vertex[ico_face[i][0]]), (ico_vertex[ico_face[i][1]]), (ico_vertex[ico_face[i][2]]) ))
    else:
        triangle_ylw = [((ico_vertex[ico_face[i][0]]), (ico_vertex[ico_face[i][1]]), (ico_vertex[ico_face[i][2]]))]
ax[0].remove()
ax[0] = fig.add_subplot(2, 2, 1, projection='3d')

tri = Poly3DCollection(triangles)
tri.set_edgecolor('k')
tri.set_facecolor(blu)
tri2 = Poly3DCollection(triangle_ylw)
tri2.set_edgecolor('k')
tri2.set_facecolor(ylw)

ax[0].add_collection(tri)
ax[0].add_collection(tri2)

lim = 1.4
ax[0].set_xlim([-lim, lim])
ax[0].set_ylim([-lim, lim])
ax[0].set_zlim([-lim, lim])
ax[0].set_box_aspect((1.0, 1.0, 1.0))

ax[0].axis('off')
#ax[0].title.set_text('(a)', y=-0.01)
ax[0].set_title('(a)', y=-0.07, fontsize=10)


### Draw and annotate the triangular face

v1 = np.array([0, 0.1])
v2 = np.array([1, 0.1])
v0 = np.array([0.5, 0.5*math.sqrt(3.0) + 0.1])

ax[1].plot([v1[0], v2[0]], [v1[1], v2[1]], 'k', linestyle="solid")
ax[1].plot([v1[0], v0[0]], [v1[1], v0[1]], 'k', linestyle="solid")
ax[1].plot([v0[0], v2[0]], [v0[1], v2[1]], 'k', linestyle="solid")

n = 4
# Triangle vertices
count = 0
v012 = []
for j in range(n+1):
    v01 = (v0*(n-j) + v1*(j)) / n
    v02 = (v0*(n-j) + v2*(j)) / n
    for k in range(j+1):
        count += 1
        if j == 0:
            v012.append(v01)
        else:
            v012.append((v01*(j-k) + v02*k) / j)

X = np.array([v1, v2, v0])
t1 = plt.Polygon(X, color=ylw)
ax[1].add_patch(t1)

neighbour = [[], [2, 4], [4], [4,7], [5,7,8], [8], [7,11], [8,11,12], [9,12,13], [13], [], [], [], [], []]
for i in range(len(v012)):
    for j in neighbour[i]:
        ax[1].plot([v012[i][0], v012[j][0]], [v012[i][1], v012[j][1]], 'b', linestyle="dotted")

for i in range(len(v012)):
    ax[1].plot(v012[i][0], v012[i][1], 'k.', markersize=8)

lim1 = -0.1
lim2 = 0.9
ax[1].set_xlim([-0.1, 1.1])
ax[1].set_ylim([-0.1, 1.1])
ax[1].axis('off')
ax[1].set_title('(b)', y=-0.07, fontsize=10)


### Plot spherical triangulation
vertex = leod.triangulation_sphere.triangulate_sphere(1.0, 4)

ax[2].remove()
ax[2] = fig.add_subplot(2, 2, 3, projection='3d')

for i in range(len(vertex)):
    for j in vertex[i].neighbour.keys():
        if j > i:
            ico_index = 0
            for k in vertex[i].neighbour[j].face:
                if k > j:
                    # Get face vertices in the appropriate format
                    x = [vertex[i].carts[0], vertex[j].carts[0], vertex[k].carts[0]]
                    y = [vertex[i].carts[1], vertex[j].carts[1], vertex[k].carts[1]]
                    z = [vertex[i].carts[2], vertex[j].carts[2], vertex[k].carts[2]]
                    verts = [list(zip(x, y, z))]
                    
                    if ylw_index in vertex[i].ico and ylw_index in vertex[j].ico and ylw_index in vertex[k].ico:
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=ylw, edgecolor='k', linewidth=0.6)
                    else:
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=blu, edgecolor='k', linewidth=0.6)
                        
                    ax[2].add_collection3d(srf)
                    ico_index += 1

limit = 0.9
ax[2].set_xlim([-limit, limit])
ax[2].set_ylim([-limit, limit])
ax[2].set_zlim([-limit, limit])
ax[2].set_box_aspect((1.0, 1.0, 1.0))
ax[2].axis('off')
ax[2].set_title('(c)', y=-0.01, fontsize=10)


### Plot ellipsoidal triangulation

ax[3].remove()
ax[3] = fig.add_subplot(2, 2, 4, projection='3d')

# Scale vertices to the surface of the ellipsoid.
for i in range(len(vertex)):
    vertex[i].carts[0] *= 3.0
    vertex[i].carts[1] *= 2.0
    vertex[i].carts[2] *= 1.0

for i in range(len(vertex)):
    for j in vertex[i].neighbour.keys():
        if j > i:
            ico_index = 0
            for k in vertex[i].neighbour[j].face:
                if k > j:
                    # Get face vertices in the appropriate format
                    x = [vertex[i].carts[0], vertex[j].carts[0], vertex[k].carts[0]]
                    y = [vertex[i].carts[1], vertex[j].carts[1], vertex[k].carts[1]]
                    z = [vertex[i].carts[2], vertex[j].carts[2], vertex[k].carts[2]]
                    verts = [list(zip(x, y, z))]
                    
                    if ylw_index in vertex[i].ico and ylw_index in vertex[j].ico and ylw_index in vertex[k].ico:
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=ylw, edgecolor='k', linewidth=0.6)
                    else:
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=blu, edgecolor='k', linewidth=0.6)
                    ax[3].add_collection3d(srf)
                    ico_index += 1
                    
limit = 1.7
ax[3].set_xlim([-limit, limit])
ax[3].set_ylim([-limit, limit])
ax[3].set_zlim([-limit, limit])
ax[3].set_box_aspect((1.0, 1.0, 1.0))
ax[3].axis('off')
ax[3].set_title('(d)', y=-0.01, fontsize=10)

fig.subplots_adjust(hspace=0.1, wspace=0.1)

if save_on:
    fig.savefig('im_icosahedral_triangulation.pdf', format='pdf', bbox_inches='tight')
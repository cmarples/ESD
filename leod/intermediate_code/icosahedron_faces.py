# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 17:07:32 2022

@author: Callum Marples

Triangulation of a sphere by subdivision of the faces of a regular icosahedron
Code based on [1]

[1] Ayad Al-Rumaithi (2022). Ellipsoid Triangulation (https://www.mathworks.com/matlabcentral/fileexchange/93505-ellipsoid-triangulation), MATLAB Central File Exchange. Retrieved September 30, 2022.
"""

import math
import numpy as np
import copy
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Sphere triangulation using a geodesic polyhedron
n = 10

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
'''
ico_vertex = [ [t, 1, 0], [-t, 1, 0], [t, -1, 0], [-t, -1, 0],
               [1, 0, t], [1, 0, -t], [-1, 0, t], [-1, 0, -t],
               [0, t, 1], [0, -t, 1], [0, t, -1], [0, -t, -1] ]
ico_vertex = np.array(ico_vertex) / den
ico_face = [ [1, 9, 5], [1, 6, 11], [3, 5, 10], [3, 12, 6], 
             [2, 7, 9], [2, 11, 8], [4, 10, 7], [4, 8, 12], 
             [1, 11, 9], [2, 9, 11], [3, 10, 12], [4, 12, 10],
             [5, 3, 1], [6, 1, 3], [7, 2, 4], [8, 4, 2],
             [9, 7, 5], [10, 5, 7], [11, 6, 8], [12, 8, 6] ]
ico_face = np.array(ico_face) - 1 # Using zero-based indexing
'''

vertex = []
face = []

# Determine triangle vertices and faces
for i in range(20):
    # Get icosahedron vertices
    v0 = ico_vertex[ico_face[i][0]]
    v1 = ico_vertex[ico_face[i][1]]
    v2 = ico_vertex[ico_face[i][2]]
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
            
    # Triangle faces
    count = 0
    face_tri = []
    for j in range(1, n+1):
        for k in range(1, j+1):
            count += 1
            no = int(j*(j-1)/2 + k) - 1
            face_tri.append([no, no+j, no+j+1])
    for j in range(2, n+1):
        for k in range(1, j):
            count += 1
            no = int(j*(j-1)/2 + k) - 1
            face_tri.append([no, no+j+1, no+1])
            
    # Merge faces and vertices
    m = len(vertex)
    for j in range(len(v012)):
        vertex.append(v012[j])
    for j in range(len(face_tri)):
        face_temp = [face_tri[j][0]+m, face_tri[j][1]+m, face_tri[j][2]+m]
        face.append(face_temp)

# Remove duplicates
q = len(vertex)
q2 = copy.copy(q)
for i in range(q-1):
    jj = i
    for j in range(i+1, q+1):
        jj += 1
        #print([i, j])
        if np.linalg.norm(vertex[i]-vertex[jj]) < 1e-5:
            del vertex[jj]
            
            #vertex[j] = np.array([0,0,0])
            q2 -= 1
            for k in range(len(face)):
                for l in range(3):
                    if face[k][l] == jj:
                        face[k][l] = i
                    elif face[k][l] > jj:
                        face[k][l] = face[k][l] - 1
            jj -= 1
        if jj >= q2-1:
            break
    if i >= q2-2:
        break

# Project onto sphere surface
for i in range(q2):
    vertex[i] /= np.linalg.norm(vertex[i])


# Get edges
edge = {}
for i in range(q2):
    edge_list = []
    for j in range(len(face)):
        if i in face[j]:
            for k in range(3):
                edge_list.append(face[j][k])
    edge[i] = []
    for j in edge_list:
        if (j not in edge[i]) and j != i:
            edge[i].append(j)

# Test triangulation by plotting
triangles = []
for i in range(len(face)):
    triangles.append(( (vertex[face[i][0]]), (vertex[face[i][1]]), (vertex[face[i][2]]) ))

fig = plt.figure(figsize=plt.figaspect(1)*2)
ax = fig.add_subplot(projection='3d', proj_type='ortho')

tri = Poly3DCollection(triangles)
tri.set_edgecolor('k')
ax.add_collection(tri)

#for i in range(q2):
#    ax.plot(vertex[i][0], vertex[i][1], vertex[i][2])

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])
ax.set_zlim([-1.1, 1.1])
plt.axis('off')

plt.show()



'''
# Plot vertices of icosahedron and unit sphere
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
# Icosahdron vertices
for i in range(12):
    ax.scatter(ico_vertex[i, 0], ico_vertex[i, 1], ico_vertex[i, 2])
# Unit sphere
th = np.linspace(0, math.pi, 25)
ph = np.linspace(0, 2*math.pi, 25)
TH, PH = np.meshgrid(th, ph)
X = 0.9 * np.sin(TH) * np.cos(PH)
Y = 0.9 * np.sin(TH) * np.sin(PH)
Z = 0.9 * np.cos(TH)
ax.plot_surface(X, Y, Z, color=[0.5, 0.5, 0.5], linewidth=0, alpha=0.1)
plt.axis('off')
'''




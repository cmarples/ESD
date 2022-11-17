# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:12:04 2022

@author: Callum Marples

Plot the triangulations used for the fast marching method, indicating the
locations of any obtuse angles.
"""

from leod.ellipsoid_shape import EllipsoidShape
from leod.triangulation_sphere import triangulate_sphere
from leod.fmm_precalculation import precalculate_grid

from math import acos
from math import pi
#import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Define shape.
shape = EllipsoidShape(3.0, 2.0, 1.0)
shape.normalise()

# Generate triangulation for the unit sphere.
n = 5  # number of triangular divisions.
vertex = triangulate_sphere(1.0, n)

# Scale vertices to the surface of the ellipsoid.
for i in range(len(vertex)):
    vertex[i].carts[0] *= shape.a_axis
    vertex[i].carts[1] *= shape.b_axis
    vertex[i].carts[2] *= shape.c_axis

# Calculate face angles.
max_angle = precalculate_grid(vertex)

# Plot faces of the triangulation
fig = plt.figure(figsize=plt.figaspect(1)*2)
ax = fig.add_subplot(projection='3d', proj_type='ortho')
rad2deg = 180.0 / pi # Conversion factor

for i in range(len(vertex)):
    for j in vertex[i].neighbour.keys():
        if j > i:
            for k in vertex[i].neighbour[j].face:
                if k > j:
                    # Find largest angle in face
                    cos_alpha_i = vertex[i].neighbour[j].face_angle[k]
                    cos_alpha_j = vertex[j].neighbour[i].face_angle[k]
                    cos_alpha_k = vertex[k].neighbour[i].face_angle[j]
                    cos_alpha = cos_alpha_i
                    if cos_alpha_j < cos_alpha:
                        cos_alpha = cos_alpha_j
                    if cos_alpha_k < cos_alpha:
                        cos_alpha = cos_alpha_k
                    alpha_max = acos(cos_alpha) * rad2deg
                    
                    # Get face vertices in the appropriate format
                    x = [vertex[i].carts[0], vertex[j].carts[0], vertex[k].carts[0]]
                    y = [vertex[i].carts[1], vertex[j].carts[1], vertex[k].carts[1]]
                    z = [vertex[i].carts[2], vertex[j].carts[2], vertex[k].carts[2]]
                    verts = [list(zip(x, y, z))]
                    
                    # Determine face colour using the angle and plot face
                    if alpha_max <= 90.0:
                        face_col = [1, 1, (alpha_max-60.0)/30.0]
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                        ax.add_collection3d(srf)
                        
                    else:
                        #face_col = [1, 1 - (alpha_max-90.0)/90.0, 1]
                        face_col = [1 - (alpha_max-90.0)/90.0, 1 - (alpha_max-90.0)/90.0, 1]
                        srf = Poly3DCollection(verts, alpha=1.0, facecolor=face_col, edgecolor='k')
                        ax.add_collection3d(srf)
                    
                    
                    




limit = shape.a_axis + 0.1
ax.set_xlim([-limit, limit])
ax.set_ylim([-limit, limit])
ax.set_zlim([-limit, limit])
plt.axis('off')

plt.show()




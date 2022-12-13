"""! @brief Defines the vertex and neighbour classes."""
##
# @ file ellipsoid_shape.py
#
# @ brief Defines the ellipsoid shape class.
#
# @ section
# - Created by Callum Marples on 04/10/2022.
# - Last modified on 11/12/2022.

import math
import numpy as np



'''
# Get the Euclidean distance between vertices i and j, updating the distance
# lists of both.
def get_distance(vertex, i, j):
    
    #if j in vertex[i].distance_to_neighbour.keys() and vertex[i].distance_to_neighbour[j] != -1.0:
    if vertex[i].neighbour[j].distance != -1.0:
        return vertex[i].neighbour[j].distance
    else:
        # Calculate distance if it has not already been found
        vertex[i].neighbour[j].distance = math.sqrt( (vertex[i].carts[0]-vertex[j].carts[0])**2.0 +
                                                     (vertex[i].carts[1]-vertex[j].carts[1])**2.0 +
                                                     (vertex[i].carts[2]-vertex[j].carts[2])**2.0 )
        vertex[j].neighbour[i].distance = vertex[i].neighbour[j].distance
        return vertex[i].neighbour[j].distance   

# Get the cosine of the angle in face (i, j, k) at face i
def get_angle(vertex, i, j, k):
    # Calculate face angle if it has not already been found
    if vertex[i].neighbour[j].face_angle[k] != 2.0:
        return vertex[i].neighbour[j].face_angle[k]
    else:
        w1 = vertex[j].carts - vertex[i].carts
        w2 = vertex[k].carts - vertex[i].carts
        a = get_distance(vertex, i, j)
        b = get_distance(vertex, i, k)
        vertex[i].neighbour[j].face_angle[k] = (w1[0]*w2[0] + w1[1]*w2[1] + w1[2]*w2[2]) / (a*b)
        vertex[i].neighbour[k].face_angle[j] = vertex[i].neighbour[j].face_angle[k]
        return vertex[i].neighbour[j].face_angle[k]
'''
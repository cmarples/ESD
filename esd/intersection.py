"""! 
@brief Ellipsoid-plane intersection routine.
@file intersection.py
@author Callum Marples
- Created by Callum Marples on 13/12/2022.
- Last modified on 30/01/2022.
"""

import math
import numpy as np

def ellipsoid_plane(shape, n_vec, q):
    """! Determine the ellipse of intersection between an ellipsoid and a plane.
         This routine implements the method given by Klein \cite Klein2012, and 
         is used in ESD to calculate triaxial taxicab distance.
    @param shape : EllipsoidShape \n
        The ellipsoid.
    @param n_vec : NumPy array of float \n
        Normal vector to the plane.
    @param q : list of float \n
        A point interior to the ellipsoid that lies within the plane.
    @return 2 element list of float \n
        The \f$a\f$ and \f$b\f$ axes of the intersecting ellipse.
    """
    
    Dq = np.array([q[0]/shape.a_axis, q[1]/shape.b_axis, q[2]/shape.c_axis])
    
    # Find in-plane vectors, r and s, such that the two are mutually orthonormal.
    if math.fabs(n_vec[0]) < 1.0e-15 and math.fabs(n_vec[1]) < 1.0e-15:
        r = np.array([1.0, 0.0, 0.0])
        s = np.array([0.0, 1.0, 0.0])
    else:
        r = np.array([n_vec[1], -n_vec[0], 0.0])
        s = np.cross(n_vec, r)
        # Normalise to one
        r /= np.linalg.norm(r)
        s /= np.linalg.norm(s)
    
    
    Dr = np.array([r[0]/shape.a_axis, r[1]/shape.b_axis, r[2]/shape.c_axis])
    Ds = np.array([s[0]/shape.a_axis, s[1]/shape.b_axis, s[2]/shape.c_axis])
    
    if np.dot(Dr, Ds) != 0.0:
        # Rotate by w to make Dr.Ds=0
        w = 0.5*math.atan2(2.0*np.dot(Dr, Ds), np.dot(Dr, Dr) - np.dot(Ds, Ds))
        cw = math.cos(w)
        sw = math.sin(w)
        r_tilde = cw*r + sw*s
        s_tilde = -sw*r + cw*s
        Dr = np.array([r_tilde[0]/shape.a_axis, r_tilde[1]/shape.b_axis, r_tilde[2]/shape.c_axis])
        Ds = np.array([s_tilde[0]/shape.a_axis, s_tilde[1]/shape.b_axis, s_tilde[2]/shape.c_axis])
    
    # Dot products
    DrDr = np.dot(Dr, Dr)
    DsDs = np.dot(Ds, Ds)
    DqDr = np.dot(Dq, Dr)
    DqDs = np.dot(Dq, Ds)
    DqDq = np.dot(Dq, Dq)
    
    # Ellipse axes
    d = DqDq - DqDr*DqDr/(DrDr) - DqDs*DqDs/(DsDs)
    a_ellipse = math.sqrt((1.0-d) / DrDr)
    b_ellipse = math.sqrt((1.0-d) / DsDs)   
    
    return [a_ellipse, b_ellipse]
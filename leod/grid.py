""" 
@brief Defines the grid class, contact distribution class and binary search routine.
@file grid.py
@author Callum Marples
- Created by Callum Marples on 12/01/2023.
- Last modified on 18/01/2022.
"""

from math import pi, ceil
from numpy import array

class Grid:
    """! The grid class.
    Contains information pertaining to a set of (polar coordinate) patches/vertices defined over 
    the surface of an ellipsoid.
    """
    def __init__(self, no_theta=19, no_phi=36):
        """! The Grid initialiser.
        @param no_theta : int (optional) \n
            The number of theta values in the grid, defaults to 19.
        @param no_phi : int (optional) \n
            The number of phi values in the grid, defaults to 36.
        """
        # Grid parameters
        self.no_theta = no_theta
        self.no_phi = no_phi
        self.no_vertices = (no_theta - 2)*no_phi + 2
        self.delta_theta = pi / (no_theta - 1)
        self.delta_phi = 2.0*pi / no_phi
        # Lists of theta and phi values
        self.theta_list = [0.0] * no_theta
        self.phi_list = [0.0] * (no_phi + 1)
        for i in range(no_theta):
            self.theta_list[i] = i * self.delta_theta
        for i in range(no_phi+1):
            self.phi_list[i] = i * self.delta_phi
            
    # Methods
    def get_vertex_index(self, theta_index, phi_index):
        """! Get scalar vertex index, given vector index (theta, phi).
        @param theta_index : int \n
            \f$\theta\f$ index of the vertex.
        @param phi_index : int \n
            \f$\phi\f$ index of the vertex.
        @return int \n
            Scalar index of the vertex
        """
        if theta_index > 0 and theta_index < self.no_theta-1:
            return 1 + phi_index + self.no_phi*(theta_index-1)
        elif theta_index == 0:
            return 0
        else: # theta_index = no_theta - 1
            return self.no_vertices - 1
        
    def get_theta_index(self, vertex_index):
        """! Get \f$\theta\f$ index, given the scalar vertex index.
        @param vertex_index : int \n
            Scalar index of the vertex.
        @return int \n
            \f$\theta\f$ index of the vertex.
        """
        if vertex_index > 0 and vertex_index < self.no_vertices-1:
            return ceil(float(vertex_index)/float(self.no_phi))
        elif vertex_index == 0:
            return 0
        elif vertex_index == self.no_vertices-1:
            return self.no_theta-1
        else:
            return self.no_theta
        
    def get_phi_index(self, vertex_index, theta_index):
        """! Get \f$\phi\f$ index, given the \f$\theta\f$ index and the scalar vertex index.
        @param vertex_index : int \n
            Scalar index of the vertex.
        @param vertex_index : int \n
            \f$\theta\f$ index of the vertex.
        @return int \n
            \f$\phi\f$ index of the vertex.
        """
        return ( vertex_index - 1 - self.no_phi*(theta_index-1) )
    
    def find_theta_index(self, th):
        """! Find the \f$\theta\f$ index closest to a given \f$\theta\f$ value,
             using the binary search method.
        @param th : float \n
            The \f$\theta\f$ value forwhich the closest index is required.
        @return int \n
            Closest \f$\theta\f$ index to the input \f$\theta\f$.
        """
        if th < pi:
            return binary_search(self.theta_list, th)
        else:
            return len(self.theta_list)-1
    
    def find_phi_index(self, ph):
        """! Find the \f$\phi\f$ index closest to a given \f$\phi\f$ value,
             using the binary search method.
        @param ph : float \n
            The \f$\phi\f$ value for which the closest index is required.
        @return int \n
            Closest \f$\phi\f$ index to the input \f$\phi\f$.
        """
        index = binary_search(self.phi_list, ph)
        if index == len(self.phi_list)-1:
            return 0
        else:
            return index
    
    def find_vertex_index(self, th, ph):
        """! Find the scalar index of the vertex closest to a given point \f$(\theta, \phi)\f$,
             using the binary search method.
         @param th : float \n
            The \f$\theta\f$ value.
        @param ph : float \n
            The \f$\phi\f$ value.
        @return int \n
            Closest vertex to the point \f$(\theta, \phi)\f$.
        """
        th_index = self.find_theta_index(th)
        ph_index = self.find_phi_index(ph)
        return self.get_vertex_index(th_index, ph_index)


def binary_search(v, x):
    """! Find index of the closest value in list v to number x
         using the binary search algorithm.
         It is assumed that the elements in v are monotonically increasing.
    @param v : list of floats \n
        List of monotonically increasing values.
    @param x : float \n
        Value for which the closest list element in v is required.
    @return int \n
        The index of the closest value in v to x.
    """
    # Initialise
    L = 0
    R = len(v)
    # Find index L such that : v[L] <= x < v[L+1]
    while L != R:
        M = ceil( 0.5*(L + R) )
        if x < v[M]:
            R = M - 1
        else:
            L = M
    # Output the index whose value is closest to x
    dL = x - v[L]
    dR = v[L+1] - x
    if dL < dR:
        return L
    else:
        return L + 1
    
    
    
class ContactDistribution():
    """! The contact distribution class.
    Contains an EllipsoidShape object, a Grid object and an array giving the 
    number of contacts per patch in the Grid. Also contains a function to add
    new contacts, given either Cartesian or polar coordinates.
    """
    def __init__(self, shape, grid):
        """! The Grid initialiser.
        @param shape : EllipsoidShape \n
            The ellipsoid.
        @param grid : Grid \n
            The grid defining the ellipsoid surface patches.
        """
        self.shape = shape
        self.grid = grid
        self.contacts = array([0] * grid.no_vertices)
    
    def add_contact(self, point):
        """! Add a new contact point to the distribution.
        @param point : list (or NumPy array) of float \n
            The contact point. Can be specified either using Cartesian \f$(x,y,z\f$
            or polar \f$(\theta, \phi\f$ coordinates.
        """
        if len(point) == 2:   # Polar coordinates
            self.update_bin(point[0], point[1])
        elif len(point) == 3: # Cartesian coordinates
            [th, ph] = self.shape.cart2polar(point[0], point[1], point[2])
            self.update_bin(th, ph)
        
    def update_bin(self, th, ph):
        """! Increments the relevant bin that contains the input polar coordinates.
        @param th : float \n
            The \f$\theta\f$ coordinate of the surface point.
        @param ph : float \n
            The \f$\phi\f$ coordinate of the surface point.
        """
        index = self.grid.find_vertex_index(th, ph)
        self.contacts[index] += 1
    
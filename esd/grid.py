"""! 
@brief Defines the grid class and a binary search routine.
@file grid.py
@author Callum Marples
- Created by Callum Marples on 12/01/2023.
- Last modified on 30/01/2022.
"""

from math import pi, ceil

class Grid:
    """! The grid class.
    In ESD, the term 'grid' refers to a set of \f$\theta\f$ and \f$\phi\f$ 
    points used to define a polar coordiante surface mesh. The grid object contains 
    the polar coordinate information in the \f$\theta\f$-\f$\phi\f$ plane. This 
    class is shape independent; it contains neither shape information, nor 
    Cartesian coordinates of surface points.
    
    Each vertex in a FmmMesh object (based on polar coordinates) can be indexed 
    by a scalar, or by a pair of integers giving the \f$\theta\f$ and \f$\phi\f$ 
    positions. Functions are defined as Grid members to convert between the two 
    and to find the index/indices of the vertex closest to a given surface point 
    (as given in polar coordinates).
    
    All angular coordiantes in this class are stored in radians.
    """
    def __init__(self, no_theta=19, no_phi=36):
        """! The Grid initialiser.
        @param no_theta : int (optional) \n
            The number of theta values in the grid, defaults to 19.
        @param no_phi : int (optional) \n
            The number of phi values in the grid, defaults to 36.
        """
        # Grid parameters
        ## Number of \f$\theta\f$ values, \f$n_\theta\f$.
        self.no_theta = no_theta
        ## Number of \f$\phi\f$ values, \f$n_\phi\f$.
        self.no_phi = no_phi
        ## Number of vertices in the grid. This is given by,
        ## \f$n_\mathrm{vertices} = (n_\theta - 2)n_\phi + 2\f$
        self.no_vertices = (no_theta - 2)*no_phi + 2
        ## Increment of \f$\theta\f$.
        self.delta_theta = pi / (no_theta - 1)
        ## Increment of \f$\phi\f$.
        self.delta_phi = 2.0*pi / no_phi
        ## List of the \f$\theta\f$ values. Ranges from 0 to \f$\pi\f$.
        self.theta_list = [0.0] * no_theta
        ## List of the \f$\phi\f$ values. Ranges from 0 to \f$2\pi\f$.
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
        @param theta_index : int \n
            \f$\theta\f$ index of the vertex.
        @return int \n
            \f$\phi\f$ index of the vertex.
        """
        return ( vertex_index - 1 - self.no_phi*(theta_index-1) )
    
    def find_theta_index(self, th):
        """! Find the \f$\theta\f$ index closest to a given \f$\theta\f$ value,
             using the binary search method.
        @param th : float \n
            The \f$\theta\f$ value, in radians, for which the closest index is required.
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
            The \f$\phi\f$ value for, in radians, which the closest index is required.
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
            The \f$\theta\f$ value, in radians.
        @param ph : float \n
            The \f$\phi\f$ value, in radians.
        @return int \n
            Closest vertex to the point \f$(\theta, \phi)\f$.
        """
        th_index = self.find_theta_index(th)
        ph_index = self.find_phi_index(ph)
        return self.get_vertex_index(th_index, ph_index)


def binary_search(v, x):
    """! Find the index of the closest value in list v to number x
         using the binary search algorithm. For further details of the algorithm, 
         see \cite Knuth1998.
         It is assumed that the elements in v are monotonically increasing.
    @param v : list of floats \n
        List of monotonically increasing values.
    @param x : float \n
        Value for which the closest list element in v is required.
    @return int \n
        The index of the closest value in v to x.
    """
    # Check extremes
    if x < v[0]:
        return 0
    elif x > v[-1]:
        return len(v)-1
    else:
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
        
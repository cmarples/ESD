# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:06:00 2022

@author: Callum Marples
"""

import numpy as np

# This class contains information for a pixel in the grid used in the fast
# marching method. An array of these objects is associated with a particular
# EllipsoidShape object.
class GeoPixel:
    def __init__(self, pixel_index, theta_index, phi_index, carts, no_pixels):
        self.pixel_index = pixel_index
        self.theta_index = theta_index
        self.phi_index = phi_index
        self.carts = np.array(carts)
        self.set_pole(no_pixels)
        self.neighbour = []
        self.neighbour_distance = []

    def set_pole(self, no_pixels):
        if self.pixel_index == 0:
            self.phi_index = 0
            self.is_north = True
            self.is_south = False
        elif self.pixel_index == no_pixels-1:
            self.phi_index = 0
            self.is_north = False
            self.is_south = True
        else:
            self.is_north = False
            self.is_south = False
    
    
    
    
    
    
    
    
    
    
    
    
    '''
    /** Default constructor. */
        GeoPixel();
        /** Constructor (input the index of the ellipsoid patch). */
        GeoPixel(int& patchNo, EllipsoidShape& E, EllipsoidPatches& P);
        /** Set pixel by specifying starting theta and phi. */
        void SetPixel(double& theta, double& phi, EllipsoidShape& E, EllipsoidPatches& P);
        /** Set pole flags based on the patch index. */
        void SetPole(int noPatches);
        /** Set Cartesian coordinates of the representative point of the pixel. */
        void SetCartesians(double& theta, double& phi, EllipsoidShape& E);
        /** Get pixel index. */
        int GetPixelIndex() { return inPixel; }

    private:
        int inPixel;                    /**< Pixel index. */
        int inTheta;                    /**< Theta index. */
        int inPhi;                      /**< Phi index. */
        std::array<double, 3> valCarts; /**< Cartesian coordinates of pixel centre. This is used to calculate Euclidean distances between pixels. */
        bool isNorth;                   /**< Is true if the pixel contains the north pole (where \f$\theta=0\f$. */
        bool isSouth;                   /**< Is true if the pixel contains the south pole (where \f$\theta=\pi\f$. */
        std::vector<int> neighbour;
        std::vector<double> neighbourDistance;
    '''
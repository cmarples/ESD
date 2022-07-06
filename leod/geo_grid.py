# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:30:04 2022

@author: Callum Marples
"""

import math
from .ellipsoid_shape import EllipsoidShape
from .geo_pixel import GeoPixel

# This class contains the set of all GeoPixel objects for a given 
# EllipsoidShape. It also includes functions for initialising the set of 
# pixels, by calculating the neighbour-to-neighbour distances.
class GeoGrid:
    # Constructor
    def __init__(self, shape=EllipsoidShape(), no_theta=19, no_phi=36):
        # Take inputs
        self.shape = shape
        self.no_theta = no_theta
        self.no_phi = no_phi
        # Calculate number of pixels and required increments
        self.no_pixels = (self.no_theta - 2)*self.no_phi + 2
        self.delta_theta = math.pi / (self.no_theta - 1)
        self.delta_phi = 2.0*math.pi / self.no_phi
        # Compute lists of theta and phi values
        self.theta_list = [0.0] * self.no_theta
        self.phi_list = [0.0] * (self.no_phi + 1)
        for i in range(self.no_theta):
            self.theta_list[i] = i * self.delta_theta
        for i in range(self.no_phi+1):
            self.phi_list[i] = i * self.delta_phi
        # Construct a list of pixel objects
        self.pixel = []
        for i in range(self.no_pixels):
            theta_index = self.find_theta_index(i)
            phi_index = self.find_phi_index(i, theta_index)
            carts = self.polars_to_cartesians( self.theta_list[theta_index],
                                               self.phi_list[phi_index] )
            self.pixel.append(GeoPixel(i, theta_index, phi_index,
                                       carts, self.no_pixels) )
            
    # Get pixel index from theta and phi indices
    def find_pixel_index(self, theta_index, phi_index):
        if theta_index > 0 and theta_index < self.no_theta-1:
            return 1 + phi_index + self.no_phi*(theta_index-1)
        elif theta_index == 0:
            return 0
        else: # theta_index = no_theta - 1
            return self.no_pixels - 1
        
    # Get theta index from pixel index
    def find_theta_index(self, pixel_index):
        if pixel_index > 0 and pixel_index < self.no_pixels-1:
            return math.ceil(float(pixel_index)/float(self.no_phi))
        elif pixel_index == 0:
            return 0
        elif pixel_index == self.no_pixels-1:
            return self.no_theta-1
        else:
            return self.no_theta
        
    # Get phi index from pixel index and theta index
    def find_phi_index(self, pixel_index, theta_index):
        return ( pixel_index - 1 - self.no_phi*(theta_index-1) )
    
    # Calculate Cartesian coordinates of a given (theta, phi) point
    # on the ellipsoid, self.shape. Return the coordinates as a list
    def polars_to_cartesians(self, theta, phi):
        sin_theta = math.sin(theta)
        return [ self.shape.a_axis * sin_theta * math.cos(phi),
                 self.shape.b_axis * sin_theta * math.sin(phi),
                 self.shape.c_axis * math.cos(theta) ]
        
    
            
    
    
'''
friend class GeoFMM;
    friend class ContactDistribution;
    public:
        GeoGrid();
        GeoGrid(EllipsoidShape& E, EllipsoidPatches& P);
        /** Get Euclidean distance between pixel i and neighbour k. If this distance has not been calculated yet, then do so. */
        double GetDistance(GeoPixel& i, int& k);
        std::vector<GeoPixel> GetPixels() { return pixel; }
        void InitialiseStickySiteParameters(double sigmaDegrees);
        double GetUpperPatchWidth() { return upperPatchWidth; }
        double GetCutoff() { return cutoffGeoDistance; }

        std::shared_ptr<EllipsoidPatches> mpPatches;    /**< Information on the set of patches. */
    private:
        std::shared_ptr<EllipsoidShape> mpShape;        /**< Ellipsoid shape on which geodesics are to be calculated. */
        std::vector<GeoPixel> pixel;                    /**< Array of pixels corresponding to the given shape. */
        double upperPatchWidth;
        double cutoffGeoDistance;
'''
"""! 
@brief Defines the ellipsoid shape class.
@file shape.py
@author Callum Marples
"""
##
# @ file ellipsoid_shape.py
#
# @ brief Defines the ellipsoid shape class.
#
# @ section
# - Created by Callum Marples on 04/07/2022.
# - Last modified on 12/12/2022.

import math
import numpy as np

class EllipsoidShape:
    """! The ellipsoid shape class.
    Contains the principal semi-axis lengths; \f$a\f$, \f$b\f$ and \f$c\f$ of an ellipsoid. 
    Also includes functions for conversion between different coordinate systems on the ellipsoid surface.
    """
    def __init__(self, a=1.0, b=1.0, c=1.0):
        """! The EllipsoidShape initialiser.
        @param a : float (optional) \n
            The \f$a\f$-axis length, defaults to 1.0.
        @param b : float (optional) \n
            The \f$b\f$-axis length, defaults to 1.0.
        @param c : float (optional) \n
            The \f$c\f$-axis length, defaults to 1.0.
        """
        self.set_axes(a, b, c)

    def is_sphere(self):
        """! Determines whether the shape is a sphere (i.e. \f$a=b=c\f$).
        @return bool \n
            True if spherical, False otherwise.
        """
        if self.a_axis == self.b_axis and self.b_axis == self.c_axis:
            return True
        else:
            return False
        
    def is_spheroid(self):
        """! Determines whether the shape is a spheroid (i.e. \f$a=b\f$).
        @return bool \n
            True if spheroidal, False otherwise.
        """
        if self.a_axis == self.b_axis and self.a_axis != self.c_axis:
            return True
        else:
            return False
        
    def set_axes(self, a, b, c, normalise=False):
        """! Set the axis lengths in the EllipsoidShape instance.
        @param a : float \n
            The \f$a\f$-axis length.
        @param b : float \n
            The \f$b\f$-axis length.
        @param c : float \n
            The \f$c\f$-axis length.
        @param normalise : bool (optional)
            If True, normalise axes using the normalise() method. Defaults to False. 
        """
        # Check that the inputs are valid
        if not (isinstance(a, (int, float)) and a > 0):
             raise ValueError(f"positive a-axis expected, got {a}")
        if not (isinstance(b, (int, float)) and b > 0):
             raise ValueError(f"positive b-axis expected, got {b}")
        if not (isinstance(c, (int, float)) and c > 0):
             raise ValueError(f"positive c-axis expected, got {c}")
        # Set axis length values
        self.a_axis = a
        self.b_axis = b
        self.c_axis = c
        if normalise == True:
            self.normalise()
    
    def normalise(self):
        """! Scale axes so that the effective radius \f$abc^{1/3}\f$ equals unity,
        or equivalently so the volume equals \f$4\pi/3\f$ (i.e. that of a unit sphere).
        This is useful for studying how some quantity varies with aspect ratio, while
        keeping the size (volume) constant.
        """
        r = (self.a_axis*self.b_axis*self.c_axis) ** (1.0/3.0)
        self.a_axis /= r
        self.b_axis /= r
        self.c_axis /= r
        
    def polar2cart(self, th, ph):
        """! Convert polar coordinates to Cartesians.
        @param th : float \n
            The \f$\theta\f$ coordinate.
        @param ph : float \n
            The \f$\phi\f$ coordinate.
        @return list \n
            The Cartesian coordinates \f$[x,y,z]\f$ of the point.
        """
        x = self.a_axis * math.sin(th) * math.cos(ph)
        y = self.b_axis * math.sin(th) * math.sin(ph)
        z = self.c_axis * math.cos(th)
        return [x, y, z]
    
    def cart2polar(self, x, y, z):
        """! Convert Cartesian coordinates to polars.
        @param x : float \n
            The \f$x\f$ coordinate.
        @param y : float \n
            The \f$y\f$ coordinate.
        @param z : float \n
            The \f$z\f$ coordinate.
        @return list \n
            The polar coordinates \f$[\theta, \phi]\f$ of the point.
        """
        th = math.acos(z / self.c_axis)
        ph = math.atan2(self.a_axis*y, self.b_axis*x)
        return [th, ph]
    
    def ellip2cart(self, be, lm):
        """! Convert ellipsoidal coordinates to Cartesians.
        @param be : float \n
            The \f$\beta\f$ coordinate.
        @param lm : float \n
            The \f$\lambda\f$ coordinate.
        @return list \n
            The Cartesian coordinates \f$[x,y,z]\f$ of the point.
        """
        a2 = self.a_axis * self.a_axis
        b2 = self.b_axis * self.b_axis
        c2 = self.c_axis * self.c_axis
        x = self.a_axis * math.cos(lm) * math.sqrt(a2 - b2*math.sin(be)**2 - c2*math.cos(be)**2) / math.sqrt(a2 - c2)
        y = self.b_axis * math.cos(be) * math.sin(lm)
        z = self.c_axis * math.sin(be) * math.sqrt(a2*math.sin(lm)**2 + b2*math.cos(lm)**2 - c2) / math.sqrt(a2 - c2)
        return [x, y, z]
    
    def cart2ellip(self, x, y, z):
        """! Convert Cartesian coordinates to ellipsoidal coordiantes.
        @param x : float \n
            The \f$x\f$ coordinate.
        @param y : float \n
            The \f$y\f$ coordinate.
        @param z : float \n
            The \f$z\f$ coordinate.
        @return list \n
            The ellipsoidal coordinates \f$[\beta, \lambda]\f$ of the point.
        """
        a2 = self.a_axis * self.a_axis
        b2 = self.b_axis * self.b_axis    # Axis lengths squared
        c2 = self.c_axis * self.c_axis
        h_x2 = a2 - c2
        h_y2 = b2 - c2    # Linear eccentricities
        h_e2 = a2 - b2
        h_x = math.sqrt(h_x2)
        
        # Initial guess (using analytical calculation)
        
        # if triaxial
        if self.a_axis > self.b_axis and self.b_axis > self.c_axis:
            
            # Construct quadratic: t^2 + k1*t + k0 = 0
            k0 = a2*b2 + a2*c2 + b2*c2 - x*x*(b2+c2) - y*y*(a2+c2) - z*z*(a2+b2)
            k1 = x*x + y*y + z*z - (a2+b2+c2)
            # Solve quadratic
            t2 = 0.5*(-k1 + math.sqrt(k1*k1 - 4*k0))
            t1 = k0 / t2
            # Ensure correct signs
            x_root = math.sqrt(a2 - t2)
            if x < 0:
                x_root *= -1
            if abs(t2 - b2) < 1e-13:
                y_root = 0
            else:
                y_root = math.sqrt(t2 - b2) 
            if y < 0:
                y_root *= -1
            if abs(t1 - c2) < 1e-13:
                z_root = 0
            else:
                z_root = math.sqrt(t1 - c2)
            if z < 0:
                z_root *= -1
            # Use roots to obtain beta and lambda
            if z == 0:
                beta = 0
            else:
                b2t1 = b2 - t1
                if abs(b2t1) < 1e-14:
                    b2t1 = 0
                beta = math.atan2(z_root, math.sqrt(b2t1))
            lmda = math.atan2(y_root, x_root)
            
        # if spheroid or sphere
        elif self.a_axis == self.b_axis or self.b_axis == self.c_axis:
            # Use simplified expressions
            beta = math.asin(z/self.c_axis)
            lmda = math.atan2(self.a_axis*y, self.b_axis*x)
            
        # Ensure lambda in correct range (-pi to pi)
        if lmda < -math.pi:
            lmda += 2*math.pi
        elif lmda > math.pi:
            lmda -= 2*math.pi
        
        x0, y0, z0 = self.ellip2cart(beta, lmda)
        d = math.sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
        # Apply numerical method
        while d > 1e-15:
            delta_r = np.array([x-x0, y-y0, z-z0])
            # Compute B and L (square roots of those in Panou 2019)
            B = math.sqrt(h_x2*math.cos(beta)**2 + h_e2*math.sin(beta)**2)
            L = math.sqrt(h_x2 - h_e2*math.cos(lmda)**2)
            # Use current guess to construct Jacobian
            J = np.array([ [-self.a_axis*h_y2*math.sin(2*beta)*math.cos(lmda) / (2*h_x*B), -self.a_axis*B*math.sin(lmda)],
                           [-self.b_axis*math.sin(beta)*math.sin(lmda), self.b_axis*math.cos(beta)*math.cos(lmda)],
                           [self.c_axis*math.cos(beta)*L/h_x, self.c_axis*h_e2*math.sin(beta)*math.sin(2*lmda) / (2*h_x*L)] ])
            corr = np.linalg.solve(J.transpose()@J, J.transpose()@delta_r) 
            beta += corr[0]
            lmda += corr[1]
            x0, y0, z0 = self.ellip2cart(beta, lmda)
            d = math.sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
        
        return [beta, lmda, d]
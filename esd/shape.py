"""!
@brief Defines the ellipsoid shape class.
@file shape.py
@author Callum Marples
- Created by Callum Marples on 04/07/2022.
- Last modified on 30/01/2022.
"""

from math import pi, sqrt, sin, cos, asin, acos, atan2
import numpy as np

class EllipsoidShape:
    """! 
    The ellipsoid shape class.
    Contains the principal semi-axis lengths; \f$a\f$, \f$b\f$ and \f$c\f$ of an 
    ellipsoid with equation,\n
    \n
    \f$\left( x/a \right)^2 + \left( y/b \right)^2 + \left( z/c \right)^2 = 1\f$,\n
    \n
    where \f$x\f$, \f$y\f$ and \f$z\f$ are Cartesian coordinates.
    Functions are included to tell whether a given shape is a sphere or spheroid, 
    based on the axis lengths. 
    Also includes functions for conversion between polar and Cartesian coordinates, 
    as well as ellipsoidal to Cartesian coordinates.
    The polar coordiantes \f$\theta\f$ and \f$\phi\f$ are defined such that,\n
    \n
    \f$ x = a \sin\theta \cos\phi \f$,\n
    \f$ y = b \sin\theta \sin\phi \f$,\n
    \f$ z = c \cos\theta \f$,
    \n    
    where \f$\theta \in [0, \pi]\f$ and \f$\phi \in [0, 2\pi)\f$. For the purposes
    of the calculations in this library, it is assumed that \f$a \geq b \geq c\f$.
    The only exception to this is a prolate spheroid, for which \f$a = b < c\f$, 
    so that the \f$c\f$-axis (from which \f$\theta\f$ is measured) is distinct.
    \n
    The ellipsoidal coordiantes \f$\beta\f$ and \f$\lambda\f$ are defined according to,\n
    \n
    \f$ x = a \cos\lambda \left( \cos^2\beta + \frac{a^2-b^2}{a^2-c^2} \sin^2\beta \right)^{1/2} \f$,\n
    \f$ y = b \cos\beta \sin\lambda\f$,\n 
    \f$ z = c \sin\beta \left( 1 - \frac{a^2-b^2}{a^2-c^2} \cos^2\lambda \right)^{1/2} \f$,\n
    \n
    where \f$\beta \in [-\pi /2, \pi /2]\f$ and \f$\lambda \in (-\pi, \pi]\f$.
    These coordinates are used to perform geodesic calculations using the 
    boundary value method of Panou \cite Panou2013.
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
        This routine will not recognise a spheroid with \f$a=c\f$ or \f$b=c\f$,
        since it is assumed here that the \f$c\f$-axis is always distinct.
        @return bool \n
            True if spheroidal, False otherwise.
        """
        if self.a_axis == self.b_axis and self.a_axis != self.c_axis:
            return True
        else:
            return False
        
    def set_axes(self, a, b, c, normalise=False):
        """! Set the axis lengths in the EllipsoidShape instance.
        This is used in the initialiser, but can also be used on an already 
        existing EllipsoidShape object.
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
        ## The \f$a\f$ semi-axis of the ellipsoid. The \f$\phi\f$ angle equals 0 or \f$\pi\f$ at this axis.
        self.a_axis = a
        ## The \f$b\f$ semi-axis of the ellipsoid. The \f$\phi\f$ angle equals \f$\pi/2\f$ or \f$3\pi/2\f$ at this axis.
        self.b_axis = b
        ## The \f$c\f$ semi-axis of the ellipsoid. The \f$\theta\f$ angle is defined from this axis.
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
            The \f$\theta\f$ coordinate, in radians.
        @param ph : float \n
            The \f$\phi\f$ coordinate, in radians.
        @return list \n
            The Cartesian coordinates \f$[x,y,z]\f$ of the point.
        """
        sin_th = sin(th)
        sin_ph = sin(ph)       
        x = self.a_axis * sin_th * sqrt(1.0 - sin_ph*sin_ph)
        if (ph > 0.5*pi and ph < 1.5*pi):
            x *= -1.0
        y = self.b_axis * sin_th * sin_ph
        z = self.c_axis * sqrt(1.0 - sin_th*sin_th)
        if (th > 0.5*pi):
            z *= -1.0
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
            The polar coordinates \f$[\theta, \phi]\f$ of the point, in radians.
        """
        th = acos(z / self.c_axis)
        ph = atan2(self.a_axis*y, self.b_axis*x)
        if ph < 0.0:
            ph += 2.0*pi
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
        x = self.a_axis * cos(lm) * sqrt(a2 - b2*sin(be)**2 - c2*cos(be)**2) / sqrt(a2 - c2)
        y = self.b_axis * cos(be) * sin(lm)
        z = self.c_axis * sin(be) * sqrt(a2*sin(lm)**2 + b2*cos(lm)**2 - c2) / sqrt(a2 - c2)
        return [x, y, z]
    
    def cart2ellip(self, x, y, z):
        """! Convert Cartesian coordinates to ellipsoidal coordiantes.
        This routine uses the numerical method given by Panou and Korakitis \cite Panou2019.
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
        h_x = sqrt(h_x2)
        
        # Initial guess (using analytical calculation)
        
        # if triaxial
        if self.a_axis > self.b_axis and self.b_axis > self.c_axis:
            
            # Construct quadratic: t^2 + k1*t + k0 = 0
            k0 = a2*b2 + a2*c2 + b2*c2 - x*x*(b2+c2) - y*y*(a2+c2) - z*z*(a2+b2)
            k1 = x*x + y*y + z*z - (a2+b2+c2)
            # Solve quadratic
            t2 = 0.5*(-k1 + sqrt(k1*k1 - 4*k0))
            t1 = k0 / t2
            # Ensure correct signs
            x_root = sqrt(a2 - t2)
            if x < 0:
                x_root *= -1
            if abs(t2 - b2) < 1e-13:
                y_root = 0
            else:
                y_root = sqrt(t2 - b2) 
            if y < 0:
                y_root *= -1
            if abs(t1 - c2) < 1e-13:
                z_root = 0
            else:
                z_root = sqrt(t1 - c2)
            if z < 0:
                z_root *= -1
            # Use roots to obtain beta and lambda
            if z == 0:
                beta = 0
            else:
                b2t1 = b2 - t1
                if abs(b2t1) < 1e-14:
                    b2t1 = 0
                beta = atan2(z_root, sqrt(b2t1))
            lmda = atan2(y_root, x_root)
            
        # if spheroid or sphere
        elif self.a_axis == self.b_axis or self.b_axis == self.c_axis:
            # Use simplified expressions
            beta = asin(z/self.c_axis)
            lmda = atan2(self.a_axis*y, self.b_axis*x)
            
        # Ensure lambda in correct range (-pi to pi)
        if lmda < -pi:
            lmda += 2*pi
        elif lmda > pi:
            lmda -= 2*pi
        
        x0, y0, z0 = self.ellip2cart(beta, lmda)
        d = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
        # Apply numerical method
        while d > 1e-15:
            delta_r = np.array([x-x0, y-y0, z-z0])
            # Compute B and L (square roots of those in Panou 2019)
            B = sqrt(h_x2*cos(beta)**2 + h_e2*sin(beta)**2)
            L = sqrt(h_x2 - h_e2*cos(lmda)**2)
            # Use current guess to construct Jacobian
            J = np.array([ [-self.a_axis*h_y2*sin(2*beta)*cos(lmda) / (2*h_x*B), -self.a_axis*B*sin(lmda)],
                           [-self.b_axis*sin(beta)*sin(lmda), self.b_axis*cos(beta)*cos(lmda)],
                           [self.c_axis*cos(beta)*L/h_x, self.c_axis*h_e2*sin(beta)*sin(2*lmda) / (2*h_x*L)] ])
            corr = np.linalg.solve(J.transpose()@J, J.transpose()@delta_r) 
            beta += corr[0]
            lmda += corr[1]
            x0, y0, z0 = self.ellip2cart(beta, lmda)
            d = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
        
        return [beta, lmda, d]
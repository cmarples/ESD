# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:02:30 2022

@author: Callum Marples
"""

import math
import numpy as np
import scipy.integrate

def bvm_dist(shape, start, end, is_radians=False, tol=1e-12, Jacobi=False, n=10000):
    
    start_temp = [0.0, 0.0]
    end_temp = [0.0, 0.0]
    # Convert to radians if input points given in degrees (assumed by default).
    if is_radians == False:
        conv = math.pi / 180.0
    else:
        conv = 1.0
    for i in range(2):
        start_temp[i] = start[i] * conv
        end_temp[i] = end[i] * conv
        
    th_0 = start_temp[0]
    ph_0 = start_temp[1]
    th_1 = end_temp[0]
    ph_1 = end_temp[1]
                   
    if th_0 == th_1 and ph_0 == ph_1:
        # Start and end points are the same => distance = 0
        if Jacobi:
            x, y, z = shape.ellip2cart(th_0, ph_0)
        else:
            x, y, z = shape.polar2cart(th_0, ph_0)
        return 0, [x, y, z]
        
    else:
        
        a2 = shape.a_axis * shape.a_axis
        b2 = shape.b_axis * shape.b_axis    # Axis lengths squared
        c2 = shape.c_axis * shape.c_axis
        h_x2 = a2 - c2
        h_y2 = b2 - c2    # Linear eccentricities
        h_e2 = a2 - b2
        
        if shape.a_axis > shape.b_axis and shape.b_axis > shape.c_axis: # Triaxial ellipsoid
            
            if Jacobi:
                beta_0, lmda_0 = th_0, ph_0
                beta_1, lmda_1 = th_1, ph_1
            else:
                # Convert to Jacobi's coordinates
                x_0, y_0, z_0 = shape.polar2cart(th_0, ph_0)
                x_1, y_1, z_1 = shape.polar2cart(th_1, ph_1)
                beta_0, lmda_0, d_0 = shape.cart2ellip(x_0, y_0, z_0)
                beta_1, lmda_1, d_1 = shape.cart2ellip(x_1, y_1, z_1)
        
            
            # First fundamental form and derivatives
            def fund_form(beta, lmda, case1):
            
                cb_sq = math.cos(beta)**2
                sb_sq = math.sin(beta)**2
                cl_sq = math.cos(lmda)**2
                sl_sq = math.sin(lmda)**2
                s2b = math.sin(2*beta)
                s2l = math.sin(2*lmda)
                
                
                num_EG = (h_y2*cb_sq + h_e2*sl_sq)
                den_B = 1 / (h_x2 - h_y2*sb_sq)
                den_L = 1 / (h_x2 - h_e2*cl_sq)
                B = (b2*sb_sq + c2*cb_sq) * den_B 
                E = num_EG * B
                L = (a2*sl_sq + b2*cl_sq) * den_L 
                G = num_EG * L
                
                B_prime = (a2*h_y2*s2b) * den_B**2
                L_prime = - (c2*h_e2*s2l) * den_L**2
                
                E_beta = B_prime * num_EG - B*h_y2*s2b
                G_beta = - L*h_y2*s2b
                E_lmda = B*h_e2*s2l
                G_lmda = L_prime * num_EG + L*h_e2*s2l
                
                E_betalmda = B_prime*h_e2*s2l
                G_betalmda = -L_prime*h_y2*s2b
                
                if case1:
                    # x system => need double beta derivatives
                    c2b = math.cos(2*beta)
                    E_2nd = 2*a2*h_y2*den_B**2*(den_B*h_y2*s2b**2 + c2b)*num_EG - 2*h_y2*(B_prime*s2b + B*c2b)               
                    G_2nd = -2*L*h_y2*c2b
                else:
                    # y system => need double lmda derivatives
                    c2l = math.cos(2*lmda)
                    E_2nd = 2*B*h_e2*c2l
                    G_2nd = 2*c2*h_e2*s2l*den_L**2*(den_L*h_e2*s2l - 1)*num_EG + 2*h_e2*(L_prime*s2l + L*c2l)
                    
                return E, E_beta, E_lmda, E_2nd, E_betalmda, G, G_beta, G_lmda, G_2nd, G_betalmda
            
            
            # Derivatives of x2 and x4
            def x_prime(x, lmda):
                
                E, E_beta, E_lmda, E_2nd, E_betalmda, G, G_beta, G_lmda, G_2nd, G_betalmda = fund_form(x[0], lmda, True)
                E2, G2 = E*E, G*G
                
                p_1 = 0.5*G_lmda/G - E_lmda/E
                p_2 = G_beta/G - 0.5*E_beta/E
                p_3 = -0.5*E_lmda/G
                
                x2_prime = p_3*x[1]**3 + p_2*x[1]**2 + p_1*x[1] + 0.5*G_beta/E
            
                # Coefficients to find x4_prime
                p_00 = 0.5*(E*G_2nd - E_beta*G_beta) / E2
                p_11 = ( 0.5*(G*G_betalmda - G_beta*G_lmda) / G2
                       - (E*E_betalmda - E_beta*E_lmda) / E2 )
                p_22 = ( (G*G_2nd - G_beta*G_beta) / G2
                       - 0.5*(E*E_2nd - E_beta*E_beta) / E2)
                p_33 = -0.5*(G*E_betalmda - E_lmda*G_beta) / G2
                # Output is x4_prime
                x4_prime = ((p_33*x[1]**3 + p_22*x[1]**2 + p_11*x[1] + p_00)*x[2]
                            + (3*p_3*x[1]**2 + 2*p_2*x[1] + p_1)*x[3])
            
                return np.array([x[1], x2_prime, x[3], x4_prime])
                
            def y_prime(y, beta):
            
                E, E_beta, E_lmda, E_2nd, E_betalmda, G, G_beta, G_lmda, G_2nd, G_betalmda = fund_form(beta, y[0], False)
                E2, G2 = E*E, G*G
                        
                q_1 = 0.5*E_beta/E - G_beta/G
                q_2 = E_lmda/E - 0.5*G_lmda/G
                q_3 = -0.5*G_beta/E
                y2_prime = q_3*y[1]**3 + q_2*y[1]**2 + q_1*y[1] + 0.5*E_lmda/G  
            
                # Coefficients to find y4_prime
                q_00 = 0.5*(E_2nd*G - E_lmda*G_lmda) / G2
                q_11 = ( 0.5*(E*E_betalmda - E_beta*E_lmda) / E2
                       - (G*G_betalmda - G_beta*G_lmda) / G2 )
                q_22 = ( (E*E_2nd - E_lmda*E_lmda) / E2
                       - 0.5*(G*G_2nd - G_lmda*G_lmda) / G2 )
                q_33 = -0.5*(E*G_betalmda - E_lmda*G_beta) / E2
                # Output is x4_prime
                y4_prime =  ( (q_33*y[1]**3 + q_22*y[1]**2 + q_11*y[1] + q_00) * y[2]
                       + (3*q_3*y[1]**2 + 2*q_2*y[1] + q_1) * y[3] )
            
                return np.array([y[1], y2_prime, y[3], y4_prime])  
            
            
            def vec_EG(beta, lmda):
                
                cb_sq = np.cos(beta)**2
                sb_sq = np.sin(beta)**2
                cl_sq = np.cos(lmda)**2
                sl_sq = np.sin(lmda)**2
                
                num_EG = (h_y2*cb_sq + h_e2*sl_sq)
                E = num_EG * (b2*sb_sq + c2*cb_sq) / (h_x2 - h_y2*sb_sq)
                G = num_EG * (a2*sl_sq + b2*cl_sq) / (h_x2 - h_e2*cl_sq)
                
                return E, G
            
            
            # Set up and solve system of ODE's
            if lmda_0 != lmda_1:
                # Integrate over lambda
                # Linear system; x1, x2, x3, x4 where:
                # 
                # x1 = beta
                # x2 = d(beta) / d(lmda)
                # x3 = d(beta) / d(beta_0_prime)
                # x4 = d(beta_prime) / d(beta_0_prime)
                # See Panou 2013
                
                h = abs(lmda_1 - lmda_0) / n
                lmda_vec = np.linspace(lmda_0, lmda_1, n+1)
                # Initial values (use sphere case)
                cb0 = math.cos(beta_0)
                cb1 = math.cos(beta_1)
                sb0 = math.sin(beta_0)
                sb1 = math.sin(beta_1)
                cw = sb0*sb1 + cb0*cb1*math.cos(lmda_1 - lmda_0)
                w = math.acos(cw)
                A = math.acos( (sb1 - cw*sb0) / (cb0*math.sin(w)) )
                
                if A != 0:
                    ddlmda = cb0 / math.tan(A)
                else:
                    ddlmda = (beta_1 - beta_0) / (lmda_1 - lmda_0)
                
                x0 = [beta_0, ddlmda, 0.0, 1.0]
                 
                x = solve_RK4(x_prime, x0, h, lmda_vec, beta_1, tol, n)
                  
                # Calculate distance along geodesic
                # Determine integrand using vectorised calculation
                beta_vec = x[:, 0]
                E, G = vec_EG(beta_vec, lmda_vec)
                I = np.sqrt(E * x[: ,1]**2 + G)
                
            else: # lmda_0 == lmda_1
                # Integrate over beta
                # Linear system; y1, y2, y3, y4 where:
                # 
                # y1 = beta
                # y2 = d(lmda) / d(beta)
                # y3 = d(lmda) / d(lmda_0_prime)
                # y4 = d(lmda_prime) / d(lmda_0_prime)
                # See Panou 2013, Section 3.2
                
                h = abs(beta_1 - beta_0) / n
                beta_vec = np.linspace(beta_0, beta_1, n+1)
                
                # Initial values
                cb0 = math.cos(beta_0)
                cb1 = math.cos(beta_1)
                sb0 = math.sin(beta_0)
                sb1 = math.sin(beta_1)
                cw = sb0*sb1 + cb0*cb1*math.cos(lmda_1 - lmda_0)
                w = math.acos(cw)
                cosA = (sb1 - cw*sb0) / (cb0*math.sin(w))
                if cosA < -1:
                    A = math.pi
                elif cosA > 1:
                    A = 0
                else:
                    A = math.acos( (sb1 - cw*sb0) / (cb0*math.sin(w)) )
                
                y0 = [lmda_0, math.tan(A) / cb0, 0.0, 1.0]
                
                y = solve_RK4(y_prime, y0, h, beta_vec, lmda_1, tol, n)
                  
                # Calculate distance along geodesic
                # Determine integrand using vectorised calculation
                lmda_vec = y[:, 0]
                E, G = vec_EG(beta_vec, lmda_vec)
                I = np.sqrt(E + G * y[:, 1]**2)
                
            #s = np.trapz(I, dx=h)
            s = scipy.integrate.simps(I, dx=h)
            
        elif (shape.a_axis == shape.b_axis and shape.b_axis != shape.c_axis) or (shape.a_axis != shape.b_axis and shape.b_axis == shape.c_axis): # Spheroid case
            
            h2 = h_x2
            
            lmda_0, lmda_1 = ph_0, ph_1
            if Jacobi:
                beta_0, beta_1 = th_0, th_1
            else:
                # Convert to Jacobi's coordinates
                beta_0 = math.asin(math.cos(th_0))
                beta_1 = math.asin(math.cos(th_1))
                
            if lmda_0 != lmda_1:
                
                # First fundamental form and derivatives
                def fund_form_spheroid(beta, lmda):
                    cb_sq = math.cos(beta)**2
                    sb_sq = math.sin(beta)**2
                    s2b = math.sin(2*beta)
            
                    c2b = math.cos(2*beta)
                    h_cos_beta_2 = h2 * cb_sq
                    
                    B = (a2*sb_sq + c2*cb_sq) / h_cos_beta_2
                    L = a2 / h2
                    
                    E = B * h_cos_beta_2
                    G = L * h_cos_beta_2
                    
                    B_prime = a2 * s2b / (h_cos_beta_2*cb_sq)
                    B_primeprime = 2*a2 * (s2b**2/cb_sq**2 + c2b/cb_sq) / h_cos_beta_2
                    
                    E_beta = B_prime*h_cos_beta_2 - B*h2*s2b
                    G_beta = -L*h2*s2b
                    
                    E_betabeta = B_primeprime*h_cos_beta_2 - 2*h2*(B_prime*s2b + B*c2b)
                    G_betabeta = -2*L*h2*c2b
                    
                    return E, E_beta, E_betabeta, G, G_beta, G_betabeta
                
                # Derivatives of system of equations
                def x_prime_spheroid(x, lmda):
                    E, E_beta, E_betabeta, G, G_beta, G_betabeta = fund_form_spheroid(x[0], lmda)
                    E2, G2 = E*E, G*G
                    
            
                    p_2 = G_beta/G - 0.5*E_beta/E
                    x2_prime = p_2*x[1]**2 + 0.5*G_beta/E
                
                    # Coefficients to find x4_prime
                    p_00 = 0.5*(E*G_betabeta - E_beta*G_beta) / E2
            
                    p_22 = ( (G*G_betabeta - G_beta*G_beta) / G2
                           - 0.5*(E*E_betabeta - E_beta*E_beta) / E2)
                    # Output is x4_prime
                    x4_prime = ( (p_22*x[1]**2 + p_00)*x[2]
                                 + 2*p_2*x[1]*x[3] )
                
                    return np.array([x[1], x2_prime, x[3], x4_prime])
                
                def vec_EG_spheroid(beta, lmda):
                    
                    cb_sq = np.cos(beta)**2
                    E = a2*np.sin(beta)**2 + c2*cb_sq
                    G = cb_sq * a2
                    
                    return E, G
                
                h = abs(lmda_1 - lmda_0) / n
                lmda_vec = np.linspace(lmda_0, lmda_1, n+1)
                # Initial values (use sphere case)
                cb0 = math.cos(beta_0)
                cb1 = math.cos(beta_1)
                sb0 = math.sin(beta_0)
                sb1 = math.sin(beta_1)
                cw = sb0*sb1 + cb0*cb1*math.cos(lmda_1 - lmda_0)
                w = math.acos(cw)
                A = math.acos( (sb1 - cw*sb0) / (cb0*math.sin(w)) )
                
                x0 = [beta_0, cb0 / math.tan(A), 0.0, 1.0]
                 
                x = solve_RK4(x_prime_spheroid, x0, h, lmda_vec, beta_1, tol, n)
                  
                # Calculate distance along geodesic
                # Determine integrand using vectorised calculation
                beta_vec = x[:, 0]
                E, G = vec_EG_spheroid(beta_vec, lmda_vec)
                I = np.sqrt(E * x[: ,1]**2 + G)
                
            else:
                # Integrate over great ellipse of constant lmda
                h = abs(beta_1 - beta_0) / n
                beta_vec = np.linspace(beta_0, beta_1, n+1)
                lmda_vec = np.array([lmda_0] * (n+1))
                I = np.sqrt(a2*np.sin(beta_vec)**2 + c2*np.cos(beta_vec)**2)
            
            #s = np.trapz(I, dx=h)
            s = scipy.integrate.simps(I, dx=h)
            
        X = shape.a_axis * np.cos(lmda_vec) * np.sqrt(abs(a2 - b2*np.sin(beta_vec)**2 - c2*np.cos(beta_vec)**2)) / np.sqrt(abs(a2 - c2))
        Y = shape.b_axis * np.cos(beta_vec) * np.sin(lmda_vec)
        Z = shape.c_axis * np.sin(beta_vec) * np.sqrt(abs(a2*np.sin(lmda_vec)**2 + b2*np.cos(lmda_vec)**2 - c2)) / np.sqrt(abs(a2 - c2))
        path_positions = np.zeros((n, 3))
        for i in range(n):
            path_positions[i] = [X[i], Y[i], Z[i]]
        
        return s, path_positions
    
def solve_RK4(fun, z0, h, t, z_end, tol, n):
    m = 0
    diff = 1
    # Solve for geodesic
    while diff > tol and m < 25:

        z = np.array((n+1) * [z0])

        for i in range(n):

            k1 = h * fun(z[i], t[i])
            k2 = h * fun(z[i] + 0.5*k1, t[i] + 0.5*h)
            k3 = h * fun(z[i] + 0.5*k2, t[i] + 0.5*h)
            k4 = h * fun(z[i] + k3, t[i] + h)
            
            z[i+1] = z[i] + (k1 + 2*k2 + 2*k3 + k4)/6.0  
            
        m += 1
        diff = z_end - z[n][0]
        # Apply correction
        z0[1] += diff / z[n][2]
        # Absolute value of difference
        diff = abs(diff)
    
    return z
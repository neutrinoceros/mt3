#-*-coding:utf-8-*-

#======================================================================
# common data and functions used by :
#   * plot_transfert_function.py
#   * fit_transfert_function.py
#======================================================================

import numpy as np

#defintion of useful functions
#======================================================================

def module(z) :
    return (z.real**2 + z.imag**2)**(1./2)

def argument(z) :
    return np.arcsin(z.imag/module(z))

#defintion of theoretical transfert function
#======================================================================

jc2d = 36525      #conversion factor from julian centuries to days
d2s  = 24*3600    #conversion factor from days to seconds

#numerical values extracted from table 1 :
#-------------------------------------------------------------------

#rigid earth flattening
eR    = 0.0032845075

#Earth's rotation rate (rad.s^-1) ---> (rad.jc^-1)
Om    = 7.292115e-5 * d2s * jc2d

#Moments of inertia (kg.m^2)
A     = 8.0115e37
Am    = 7.0999e37
Af    = 9.0583e37

#ellipticities
e     = 3.257e-3

#Compliances
kappa = 1.039e-3
gamma = 1.965e-3


#other numerical values 
#-------------------------------------------------------------------
sigCW   =  1.816829e-7 * d2s * jc2d #(rad.s^-1) ---> (rad.jc^-1)
sigNDFW = -7.308004e-5 * d2s * jc2d #(rad.s^-1) ---> (rad.jc^-1)



#theoretical transfert function
#--------------------------------------------------
def th_T(sig0,p1=np.zeros(4)) :
    """defined as eq 54 in "Drilling to the center of the Earth with VLBI" """
    sig = sig0 - Om
    res = - (sig-eR*Om)/(eR*Om) * (
                                   ((kappa+p1[0]) - Af/A * gamma) 
                                   - Am/A * (sigCW+p1[1])/(sig - (sigCW+p1[1]))
                                   + Af/A * ((e+p1[2]) - gamma) * ((sigNDFW+p1[3]) + Om)/(sig - (sigNDFW+p1[3]))
                                   )
    return res


#partial derivatives of theoretical transfert function
#------------------------------------------------------

def dTkappa(sig0) :
    sig = sig0 - Om
    res = - (sig-eR*Om)/(eR*Om)
    return res

def dTsigCW(sig0) :
    sig = sig0 - Om
    res = (sig-eR*Om)/(eR*Om) * Am/A * sig/(sig - sigCW)**2
    return res

def dTe(sig0) :
    sig = sig0 - Om
    res = - (sig-eR*Om)/(eR*Om) * Af/A * (sigNDFW + Om)/(sig - sigNDFW)
    return res

def dTsigNDFW(sig0) :
    sig = sig0 - Om
    res = - (sig-eR*Om)/(eR*Om) * Af/A * (e - gamma) * (sig)/(sig - sigNDFW)**2
    return res


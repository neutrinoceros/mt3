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
ef    = 2.547e-3

#Compliances
kappa = 1.039e-3
gamma = 1.965e-3
beta  = 6.160e-4

#other numerical values 
#-------------------------------------------------------------------
sigCW   =  1.816829e-7 * d2s * jc2d #(rad.s^-1) ---> (rad.jc^-1)
sigNDFW = -7.308004e-5 * d2s * jc2d #(rad.s^-1) ---> (rad.jc^-1)

efb     = ef - beta #home made variable

#theoretical transfert function
#--------------------------------------------------
def th_T(sig0,p1=np.zeros(4)) :
    """defined as eq 54 in "Drilling to the center of the Earth with VLBI" """
    sig = sig0 - Om
    kkappa = kappa + p1[0]
    ggamma = gamma + p1[1]
    ee     = e     + p1[2]
    eefb   = efb   + p1[3]
    res = - (sig-eR*Om)/(eR*Om) * (
                                   kkappa - Af/A * ggamma
                                   - Om*(ee - kkappa) / (sig - A/Am *Om*(ee - kkappa))
                                   - Af/Am * Om * (ee - ggamma) * eefb / (sig + Om *(1. + A/Am * eefb))
    )
    return res


#partial derivatives of theoretical transfert function
#------------------------------------------------------

def dTkappa(sig0) :
    sig = sig0 - Om
    res = (eR*Om-sig)/(eR*Om) * (1. 
                                 + Om/(sig-Om*(e-kappa))
                                 + Om**2*(e-kappa)/(sig-Om*(e-kappa))**2
    )
    return res

def dTgamma(sig0) :
    sig = sig0 - Om
    res = (eR*Om-sig)/(eR*Om) * ( - Af/A 
                                  + Af/Am * Om*efb/(sig + Om*(1. + A/Am*efb))
    )
    return res

def dTe(sig0) :
    sig = sig0 - Om
    res = (eR*Om - sig)/(eR*Om) * ( - Om/(sig-Om*(e-kappa)) 
                                    - Om**2*(e-kappa)/(sig-Om*(e-kappa))**2 
                                    - Af/Am * Om*efb/(sig+Om*(1. + A/Am*efb))
    )
    return res

def dTefb(sig0) :
    sig = sig0 - Om
    res = (eR*Om-sig)/(eR*Om) * ( -Af/Am * Om*(e-gamma)/(sig+Om*(1. + A/Am *efb))
                                  + A*Af/Am**2 * Om**2*(e-gamma)*efb / (sig + Om *(1. + A/Am*efb))**2
    )
    return res

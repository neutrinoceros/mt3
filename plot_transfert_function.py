#-*-coding:utf-8-*-

#======================================================================
# script used to plot a transfert function etaC/etaR
# where :
#       * etaC is complex nutation amplitudes (MHB) got from 
#         ondes.txt + corrections from amplitudes.dat (fitted parameters)
#       * etaR is complex nutation amplitudes from REN (Rigid Earth)
#         columns 17 and 18 in ondes.txt
#
# NOTE : th_T does not show a resonance where it is supposed to ...
#         fix it plz
#
#======================================================================

import numpy as np
import pylab as pl

#defintion of usefull functions
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


def th_T(omega) : #fixme
    """defined as eq 54 in "Drilling to the center of the Earth with VLBI" """
    sigprime = omega/(2.0*np.pi)
    sig = sigprime - Om
    res = - (sig-eR*Om)/(eR*Om) * (kappa - Af/A * gamma - Am/A * sigCW/(sig - sigCW) + Af/A * (e - gamma) * (sigNDFW + Om)/(sig - sigNDFW))
    return res


#data loading
#======================================================================
ondes_tab = np.loadtxt("data/ondes.txt",comments='%',usecols=range(16,20))
corrs_tab = np.loadtxt("amplitude.dat")

etaR = ondes_tab[:,0] + 1j * ondes_tab[:,0]
etaC = (ondes_tab[:,2] + corrs_tab[:,0]) + 1j * ( ondes_tab[:,3] + corrs_tab[:,1] )
transfert = np.array(etaC/etaR)

mod = module(transfert)
phi = argument(transfert)
oms = corrs_tab[:,2]


#def of theoretical values 
#----------------------------------------
oms_th = np.arange(min(oms),max(oms),1e-2)
transfert_th = th_T(oms_th)
mod_th       = module(transfert_th)
phi_th       = argument(transfert_th)


#plot script
#======================================================================

fig,axes = pl.subplots(nrows=2)

axes[1].set_xlabel(r"$\omega$",size=20)
axes[1].set_ylabel(r"$\phi$"  ,size=20)
axes[0].set_ylabel(r"$r$"     ,size=20)
axes[0].set_ylim(min(mod)*0.9,max(mod)*1.1)

pl.ion()
pl.show()

#calculated tranfert function
#----------------------------------------
axes[0].scatter(oms,mod   ,marker="+")
axes[1].scatter(oms,phi,marker='+',color='r')

#theoretical tranfert function
#----------------------------------------
axes[0].plot(oms_th,mod_th,color='m')
#axes[1].scatter(oms,phi_th,marker="*",color='m')   #useless : th_T returns reals, not complex

#theoretical resonances
#----------------------------------------
axes[0].plot(sigCW/(2*np.pi)*np.ones(2),[min(mod)/5,max(mod)*2], ls='--',c='k',lw=2)

pl.draw()
pl.ioff()

raw_input()

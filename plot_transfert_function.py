#-*-coding:utf-8-*-

#======================================================================
# script used to plot a transfert function etaC/etaR
# where :
#       * etaC is complex nutation amplitudes (MHB) got from 
#         ondes.txt + corrections from amplitudes.dat (fitted parameters)
#       * etaR is complex nutation amplitudes from REN (Rigid Earth)
#         columns 17 and 18 in ondes.txt
#
#
# also fits some parameters of the transert function :
#  * kappa (complex)
#  * sigCW and sigNDFW (complex)
#  * e (real)
#
# Note : e should be real but is treated as a complex for now (much simpler)
#======================================================================

import numpy as np
import pylab as pl
import numpy.linalg as la

from transfert import *

LW    = 1.5
ALPHA = 0.7
#data loading
#======================================================================
ondes_tab = np.loadtxt("data/ondes.txt",comments='%',usecols=range(16,20))
corrs_tab = np.loadtxt("amplitude.dat")

etaR = ondes_tab[:,0] + 1j * ondes_tab[:,0]
etaC = (ondes_tab[:,2] + corrs_tab[:,0]) + 1j * ( ondes_tab[:,3] + corrs_tab[:,1] )
transfert = np.array(etaC/etaR)

mod  = module(transfert)
phi  = argument(transfert)
sigs = corrs_tab[:,2]


#resol obs = M*p1 (find p1)
#======================================================================

obs    = transfert - th_T(sigs)
M      = np.zeros((len(obs),4))
M[:,0] = dTkappa(sigs)
M[:,1] = dTgamma(sigs)
M[:,2] = dTe(sigs)
M[:,3] = dTefb(sigs)

p1 = la.lstsq(M,obs)[0]
print p1

#def of theoretical values 
#----------------------------------------
sigMIN       = -4e4
sigMAX       = -sigMIN 
sigs_th      = np.arange(sigMIN,sigMAX,1e-1)
transfert_th = th_T(sigs_th)
mod_th       = module(transfert_th)
phi_th       = argument(transfert_th)


#def of theoretical values  (postfit)
#----------------------------------------
transfert_pf = th_T(sigs_th,p1)
mod_pf       = module(transfert_pf)
phi_pf       = argument(transfert_pf)


#plot script
#======================================================================


fig,axes = pl.subplots(nrows=2,sharex=True,squeeze=True)
pl.subplots_adjust(hspace=.001,bottom=.1)
for ax in axes :
    ax.xaxis.set_tick_params(width=1.3)
    ax.yaxis.set_tick_params(width=1.3)

axes[0].set_xlim(sigMIN,sigMAX)
axes[1].set_xlabel(r"$\sigma$",size=20)
axes[1].set_ylabel(r"arg$(T)$",size=20)
axes[0].set_ylabel(r"$|T|$"   ,size=20)
axes[0].set_ylim(0,2)

pl.ion()
pl.show()

#calculated tranfert function
#----------------------------------------
axes[0].scatter(sigs,mod,marker="+",s=50,lw=1.3)
axes[1].scatter(sigs,phi,marker="+",s=50,lw=1.3)

#theoretical tranfert function
#----------------------------------------
axes[0].plot(sigs_th,mod_th,color='m',lw=LW,alpha=ALPHA,label='a priori')
axes[0].plot(sigs_th,mod_pf,color='r',lw=LW,alpha=ALPHA,label='post-fit')
#axes[1].plot(sigs_th,phi_th,marker="*",color='m')   #useless : th_T returns reals, not complexs
axes[1].plot(sigs_th,phi_pf,color='r',lw=LW,alpha=ALPHA)   

#theoretical asymptotes
#----------------------------------------
#axes[0].plot((sigNDFW+Om)*np.ones(2),[min(mod)/5,max(mod)*2], ls='--',c='k')
#axes[0].plot(sigNDFW*np.ones(2),[min(mod)/5,max(mod)*2], ls='--',c='k',lw=2)

#xxx = np.arange(min(sigs),max(sigs))
#axes[0].plot(xxx,(Am*sigCW/(A*eR*Om))*np.ones(len(xxx)), ls='--',c='k')

axes[0].legend()
pl.draw()
pl.ioff()

raw_input()

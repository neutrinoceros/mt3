#-*-coding:utf-8-*-

#======================================================================    
# script used to plot a transfert function etaC/etaR
# where :
#       * etaC is complex nutation amplitudes (MHB) got from 
#         ondes.txt + corrections from amplitudes.dat (fitted parameters)
#       * etaR is complex nutation amplitudes from REN (Rigid Earth)
#         columns 17 and 18 in ondes.txt
#
#======================================================================    


import numpy as np
import pylab as pl

#data loading
#======================================================================    
ondes_tab = np.loadtxt("data/ondes.txt",comments='%',usecols=range(16,20))
corrs_tab = np.loadtxt("amplitude.dat")

# etaR = np.zeros((42,2))
# etaC = np.zeros((42,2))

# etaR[:,0] = ondes_tab[:,0]
# etaR[:,1] = ondes_tab[:,1]
# etaC[:,0] = ondes_tab[:,2] + corrs_tab[:,0]
# etaC[:,1] = ondes_tab[:,3] + corrs_tab[:,1]
#Tx = etaC[:,0]/etaR[:,0]
#Ty = etaC[:,1]/etaR[:,1]

etaR = ondes_tab[:,0] + 1j * ondes_tab[:,0]
etaC = (ondes_tab[:,2] + corrs_tab[:,0]) + 1j * ( ondes_tab[:,3] + corrs_tab[:,1] )
transfert = np.array(etaC/etaR)

mod = (transfert.real**2 + transfert.imag**2)**(1./2)
phi = np.arcsin(transfert.imag/mod)
amp = np.arange(len(transfert)) #TMP

#plot script
#======================================================================    

fig,axes = pl.subplots(nrows=2)
pl.ion()
pl.show()

axes[0].plot(amp,mod,marker="+")
axes[1].plot(amp,phi,marker='+',color='r')
pl.draw()
pl.ioff()

raw_input()

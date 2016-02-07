import numpy as np
import scipy as sc
import scipy.fftpack as fftp
import matplotlib.pyplot as plt
import interpol as intp
import fft_processing as ft

#===================================================================================================#
#script to interpolate with the module interpol
#then we apply the fft to the interpolated data in order to fine the resonance frequence of nutation
#===================================================================================================#


data=np.genfromtxt("../data/opa2015a.eops",usecols=(0,4,5,9,10))
data[:,0]=data[:,0]-data[0,0]
#data contain the result of VLBI session
#data[:,0] : time in julian days
#data[:,1] : Celestial pole offset dX wrt IAU 2006
#data[:,2] : Celestial pole offset dY wrt IAU 2006
#data[:,3] : Formal uncertainty of celestial pole offset dX
#data[:,4] : Formal uncertainty of celestial pole offset dY


#trying to kill some aberant point
dX_array=data[:,1] 
dY_array=data[:,2]
dX_array=dX_array-np.average(dX_array)
dY_array=dY_array-np.average(dY_array)

dX_sigma=data[:,3]
dY_sigma=data[:,4]

data_time=data[:,0]

#sample step
sample_step=15 #days

#interpolation of the data
interpolated_time_X,interpolated_dX=intp.interpolation_function(data_time,dX_array,dX_sigma,interp_time_step=sample_step)
interpolated_time_Y,interpolated_dY=intp.interpolation_function(data_time,dY_array,dY_sigma,interp_time_step=sample_step)

print interpolated_time_X[-1],interpolated_dX[-1]
raw_input()

#temporal offset put to zero
interpolated_time_X=interpolated_time_X[:]-interpolated_time_X[0]
interpolated_time_Y=interpolated_time_Y[:]-interpolated_time_Y[0]

print "nbr temps interp X =",len(interpolated_time_X)
print "nbr valeur signal X =",len(interpolated_dX)
print "nbr temps interp Y =",len(interpolated_time_Y)
print "nbr valeur signal Y =",len(interpolated_dY)


#fft processing
frequence_X,spectral_X=ft.processing_fft(interpolated_time_X,interpolated_dX,sample_step)
frequence_Y,spectral_Y=ft.processing_fft(interpolated_time_Y,interpolated_dY,sample_step)

# final time to interpolate the time
# final_time=(int(data[-1,0]/sample_step)+1)*sample_step
# print "data end",data[-1,0]
# print "final_time",final_time

#interpolation of data
# interpolated_time=np.arange(0.,final_time+sample_step,float(sample_step))
# interpolated_dX=np.interp(interpolated_time,data[:,0],dX_array)
# interpolated_dY=np.interp(interpolated_time,data[:,0],dY_array)
# Nbr_of_time=np.size(interpolated_time)
# print "end of interpolated time",interpolated_time[-1]
# print "nbr of point", Nbr_of_time

# # complex fast fourier transform normalized
# Cplx_fft_dX=fftp.fft(interpolated_dX)
# Cplx_fft_dY=fftp.fft(interpolated_dY)

# fft module processing
# N.B. the first part of tab is the conjugate of the second part
# Mod_fft_dX=np.abs(Cplx_fft_dX[:Nbr_of_time/2])*2/Nbr_of_time
# Mod_fft_dY=np.abs(Cplx_fft_dY[:Nbr_of_time/2])*2/Nbr_of_time

# frequential axis
# axe_f=np.linspace(0.,1./(2.*sample_step),Nbr_of_time/2)

print "freq",len(frequence_X)
print "spectral",len(spectral_X)

# # ploting result
plt.plot(frequence_X,spectral_X,label='nutation dX')
# plt.plot(1./axe_f,Mod_fft_dY,label='nutation dY')
plt.xlabel('jour')
plt.ylabel('puissance normalise')
plt.legend(loc='best')

#checking plot
# plt.plot(data[:,0],data[:,1],'b')
# plt.plot(interpolated_time,interpolated_dX,label='dX')
# plt.plot(interpolated_time,interpolated_dY,label='dY')
# plt.legend(loc='best')


plt.show()

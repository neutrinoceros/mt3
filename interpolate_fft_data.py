import numpy as np
import scipy as sc
import scipy.fftpack as fftp
import matplotlib.pyplot as plt

data=np.genfromtxt("data/opa2015a.eops",usecols=(0,4,5,9,10))
data[:,0]=data[:,0]-data[0,0]
#cata contain the result of VLBI session
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


#sample step
sample_step=7 #days

#final time to interpolate the time
final_time=(int(data[-1,0]/sample_step)+1)*sample_step
print "data end",data[-1,0]
print "final_time",final_time

#interpolation of data
interpolated_time=np.arange(0.,final_time+sample_step,float(sample_step))
interpolated_dX=np.interp(interpolated_time,data[:,0],dX_array)
interpolated_dY=np.interp(interpolated_time,data[:,0],dY_array)
Nbr_of_time=np.size(interpolated_time)
print "end of interpolated time",interpolated_time[-1]
print "nbr of point", Nbr_of_time

# complex fast fourier transform normalized
Cplx_fft_dX=fftp.fft(interpolated_dX)
Cplx_fft_dY=fftp.fft(interpolated_dY)

# fft module processing
# N.B. the first part of tab is the conjugate of the second part
Mod_fft_dX=np.abs(Cplx_fft_dX[:Nbr_of_time/2])*2/Nbr_of_time
Mod_fft_dY=np.abs(Cplx_fft_dY[:Nbr_of_time/2])*2/Nbr_of_time

# frequential axis
axe_f=fftp.fftfreq(len(interpolated_time),sample_step)


# # ploting result
plt.plot(axe_f[0:len(Mod_fft_dX)],Mod_fft_dX,label='nutation dX')
plt.plot(axe_f[0:len(Mod_fft_dX)],Mod_fft_dY,label='nutation dY')
plt.xlabel('frequence 1/jour')
plt.ylabel('puissance normalise')
# plt.legend(loc='best')

#checking plot
# plt.plot(data[:,0],data[:,1],'b')
# plt.plot(interpolated_time,interpolated_dX,label='dX')
# plt.plot(interpolated_time,interpolated_dY,label='dY')
plt.legend(loc='best')


plt.show()

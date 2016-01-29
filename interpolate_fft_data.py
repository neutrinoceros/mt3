import numpy as np
import scipy as sc
import scipy.fftpack as fftp
import matplotlib.pyplot as plt

data =np.genfromtxt("data/opa2015a.eops",usecols=(0,4,5,9,10))
data[:,0]=data[:,0]-data[0,0]
#cata contain the result of VLBI session
#data[:,0] : time in julian days
#data[:,1] : Celestial pole offset dX wrt IAU 2006
#data[:,2] : Celestial pole offset dY wrt IAU 2006
#data[:,3] : Formal uncertainty of celestial pole offset dX
#data[:,4] : Formal uncertainty of celestial pole offset dY

#sample step
sample_step=3 #days

#final time to interpolate the time
final_time=(int(data[-1,0]/sample_step)+1)*sample_step
# print "derniertemps data",data[-1,0]
# print "final time",final_time

#interpolation of data
interpolated_time=np.arange(0,final_time+sample_step,float(sample_step))
print "dernier temps interpoller",interpolated_time[-1] 
print "longeur tab intrerpoller",len(interpolated_time)
print "temps final / sample",final_time*1./sample_step
print "tab temp interpole", interpolated_time
interpolated_dX=np.interp(interpolated_time,data[:,0],data[:,1])
interpolated_dY=np.interp(interpolated_time,data[:,0],data[:,2])


#checking plot
# plt.plot(data[:,0],data[:,1],'b')
# plt.plot(interpolated_time,interpolated_dX,'r x')
# plt.show()


# complex fast fourier transform normalized
Cplx_fft_dX=fftp.fft(interpolated_dX)/np.size(interpolated_dX)
Cplx_fft_dY=fftp.fft(interpolated_dY)/np.size(interpolated_dY)

# fft module processing
# N.B. the first part of tab is the conjugate of the second part
Mod_fft_dX=np.abs(Cplx_fft_dX[0:len(Cplx_fft_dX)/2])*2
Mod_fft_dY=np.abs(Cplx_fft_dY[0:len(Cplx_fft_dY)/2])*2

# frequential axis
axe_f=fftp.fftfreq(len(interpolated_time),sample_step)


# ploting result
plt.plot(axe_f[0:len(Mod_fft_dX)],Mod_fft_dX,label='nutation dX')
plt.plot(axe_f[0:len(Mod_fft_dX)],Mod_fft_dY,label='nutation dY')
plt.xlabel('frequence 1/jour')
plt.ylabel('puissance normalise')
plt.legend(loc='best')
plt.show()

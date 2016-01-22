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
print "derniertemps",data[-1,0]
print "final time",final_time

#interpolation of data
interpolated_time=np.arange(0,final_time+sample_step,float(sample_step))
print "dernier temps interpoller",interpolated_time[-1] 
print "longeur tab intrerpoller",len(interpolated_time)
print "temps final * sample",final_time*sample_step
interpolated_dX=np.interp(interpolated_time,data[:,0],data[:,1])
interpolated_dY=np.interp(interpolated_time,data[:,0],data[:,2])


#checking plot
plt.plot(data[:,0],data[:,1],'b')
plt.plot(interpolated_time,interpolated_dX,'r x')
# plt.show()

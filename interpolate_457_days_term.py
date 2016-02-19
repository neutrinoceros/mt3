#-*-coding:utf-8-*-
#====================================================================================================
# interpolation script to get interpolate the amplitude of the 457 days term calculated in 
# 457_days_ampl.dat at VLBI observations dates in data/opa2015a.eops, in order to kill it later on
#
# NB : this is an interpolation-only script, it does not do any extrapolation, so you might expect
# issues if you try to evaluate the amplitude for early and late dates as we don't have calculated
# values for it near the borders of the time interval covered by VLBI sessions.
# If you apply this script to the whole data nonetheless, you will get pathological (flag is -1000)
# values where extrapolation would be recquiered. We can simply ignore those margins and work with a 
# slightly narrower data sample.
#====================================================================================================

import numpy as np
# import pylab as pl
# import matplotlib.pyplot as plt


#interpolation function
#====================================================================================================
def interp_func(t,dim) :
    if t<time[0] or t>time[-1] :
        print "error, t must be in a precise interval"
        return -1000
    else :
        if dim==1 :
            amp = xamp
        elif dim==2 :
            amp = yamp
        else :
            print "error, dim must be 1 or 2"
            return -1000
        i = 0
        while(time[i]<t) :
            i+=1
        dates = [time[i-1],time[i]]
        amps  = [amp[i-1],  amp[i]]
        p = np.poly1d(np.polyfit(dates,amps,1))
        return p(t)


#application and data saving
#====================================================================================================
# xamp_inter = np.array([interp_func(t,dim=1) for t in time_d])
# yamp_inter = np.array([interp_func(t,dim=2) for t in time_d])

# tab_final = np.column_stack([time_d,xamp_inter,yamp_inter])
# np.savetxt('interpolated_457d_amps.dat',tab_final,delimiter='    ')

#data loading
#====================================================================================================
tab   = np.loadtxt('457_days_ampl.dat')
time,xamp,yamp = tab[:,0],tab[:,1],tab[:,2]

time_d = np.loadtxt('data/opa2015a.eops',usecols=[0])



#trying with the interne function 
#================================
x_interp=np.interp(time_d,time,xamp)
y_interp=np.interp(time_d,time,yamp)

#saving value
#============
tab_final = np.column_stack([x_interp,y_interp])
np.savetxt('interpolated_457d_amps.dat',tab_final,delimiter='    ')




# plt.plot(time,xamp,'r-',label='data')
# plt.scatter(time_d,x_interp,label='interp',alpha=0.1)
# plt.legend(loc='best')
# plt.show()




#for test/visualization only
#====================================================================================================
# time_i = np.arange(min(time),max(time),1e-1)
# xamp_i = np.array([interp_func(t) for t in time_i])

# fig,ax = pl.subplots()
# pl.ion()
# pl.show()
# ax.scatter(time,xamp)
# # ax.scatter(time_d,time_d)
# ax.plot(time_i,xamp_i)
# pl.draw()
# pl.ioff()
# raw_input()

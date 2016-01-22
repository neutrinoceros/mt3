#-*-coding:utf-8-*-

#======================================================================    
# script used to interpolate by different methods our signal
# in order to obtain a regular time-step version of it
# 
# work in progress
# to do :
#     * numerically estimate which method works best
#     * export functions to a python module
#======================================================================    

import numpy as np
import pylab as pl


#======================================================================    
#                           functions def
#======================================================================    


def set_zeros2ones(w) :
    #used to avoid infinite weights that would occur when error bar is zero
    for i in range(len(w)) :
        if w[i] == 0 :
            w[i] = 1.
    w = w**(-1)


def comp_mean(xdata,ydata,ysigma,tmin,h=7.0) :
    #computes mean y-value of a ydata set within the interval [tmin, tmin+h] (in days)
    #xdata are times at which y has been measured
    #ysigma are standard errors over ydata used to weight the points

    #this method assumes mean value to be that of a linear function fited over 
    #our data points in the delimimited interval
    w = np.array(ysigma)**2
    set_zeros2ones(w)
    a,b = np.polyfit(xdata,ydata,1,w=w)
    ymean = a/2 * (h + 2*tmin) + b
    return ymean

def ponderateur(t,N,o,p) :
    #arbitrary weight function
    return N*np.sin(t*o+p)


def comp_mean2(xdata,ydata,ysigma,tmin,h=7.0) :
    #this second method computes mean value at the middle of the interval to be closer to y values 
    #known near this date than to y values known at the interval's boundaries

    x    = np.array(xdata)
    y    = np.array(ydata)
    om   = np.pi/h
    phi  = tmin*om
    Norm = np.pi/h *(np.cos(om*tmin+phi) - np.cos(om*(tmin+h)+phi))
    w = np.array(ysigma)
    set_zeros2ones(w)
    w = w**-2
    w *= ponderateur(x,Norm,om,phi)
    
    ymean = np.sum(y*w)/np.sum(w)
    return ymean


def get_mean_signal(time,xpol,ypol,sigxpol,sigypol,step=7.,ignore=200,method=1) :
    #computes mean signal for x and y at the same time according to one method or the other
    #ignore is used to avoid difficulties encountered for oldest data points when low sampling rates 
    #did not allow avering over short periods of time

    if   method == 1 :
        cm = comp_mean
    elif method == 2 :
        cm = comp_mean2

    tmean = []
    xmean = []
    ymean = []
    tmin = time[ignore]
    i = 0
    while tmin+step <= time[-1] :
        tdata,xdata,xsigma,ydata,ysigma = [],[],[],[],[]
        while time[i] < tmin + step :
            tdata.append(time[i])
            xdata.append(xpol[i])
            xsigma.append(sigxpol[i])
            ydata.append(ypol[i])
            ysigma.append(sigypol[i])
            i += 1
        tm = tmin + step/2.0
        xm = cm(tdata,xdata,xsigma,tmin,step)
        ym = cm(tdata,ydata,ysigma,tmin,step)
        tmean.append(tm)
        xmean.append(xm)
        ymean.append(ym)
        tmin += step
    return tmean,xmean,ymean


def plot_signalVSmean(ax,t,xdata,sigx,tm,xm) :
    #simple plotting function
    step = tmean[1] - tmean[0]
#    xticks = tm
#    ax.set_xticks(xticks)
    ax.errorbar(time,xdata,yerr=sigx,fmt='+',alpha=0.6,markersize=1)
#    ax.scatter(t,xdata,marker='o',s=1.,alpha=0.6)
    ax.scatter(tm,xm,s=100,lw=2,marker='+',color='r',alpha=0.8)


#======================================================================    
#                           plot script
#======================================================================    


#data loading
filename = 'data/opa2015a.eops'
tab      = np.genfromtxt(filename,usecols=(0,4,5,9,10))

time,xpol,ypol,sigxpol,sigypol = tab[:,0],tab[:,1],tab[:,2],tab[:,3],tab[:,4]

#data processing
tmean,xmean,ymean = get_mean_signal(time,xpol,ypol,sigxpol,sigypol,step=20.,ignore=200,method=2)


#plotting
pl.ion()
fig,axes =  pl.subplots(nrows=2)

pl.show()
plot_signalVSmean(axes[0],time,xpol,sigxpol,tmean,xmean)
plot_signalVSmean(axes[1],time,ypol,sigypol,tmean,ymean)


pl.draw()
pl.ioff()

raw_input()

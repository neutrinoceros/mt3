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
import numpy.ma as ma
import pylab as pl
import math

#======================================================================    
#                           functions def
#======================================================================    

def get_mask(signal,err) :
    #finds invalid values in the signal and return the corresponding mask

    mask = []
    for i in range(len(signal)) :
        if signal[i] == 0.0 or err[i] == 0.0 :
            mask.append(1)
        else :
            mask.append(0)
    return mask


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


def ponderateur(t,om) :
    #arbitrary weight function
    return (1. + np.sin(t*om+3*np.pi/2))


def comp_mean2(time,xdata,xsigma,tmin,h=7.0) :
    #this second method computes mean value at the middle of the interval to be closer to y values 
    #known near this date than to y values known at the interval's boundaries

    t    = np.array(time)
    x    = np.array(xdata)
    # om   = np.pi/h
    # phi  = tmin*om
    # Norm = np.pi/h *(np.cos(om*tmin+phi) - np.cos(om*(tmin+h)+phi))

    om   = 2*np.pi/h
    
    w = np.array(xsigma)
    w = w**-2
#    w *= ponderateur(x,Norm,om,phi)
#    w *= ponderateur(x,om)

    xmean = np.sum(x*w)/np.sum(w)
    return xmean


def get_mean_signal(time,xpol,sigxpol,mask,step=20.,ignore=200,method=2) :
    #computes mean signal for x and y at the same time according to one method or the other
    #ignore is used to avoid difficulties encountered for oldest data points when low sampling rates 
    #did not allow avering over short periods of time

    if   method == 1 :
        cm = comp_mean
    elif method == 2 :
        cm = comp_mean2

    mtime     = ma.masked_array(time,mask)
    mxpol     = ma.masked_array(xpol,mask)
    msigxpol  = ma.masked_array(sigxpol,mask)
    
    #only take non masked values
    mtime    = mtime[~mtime.mask]
    mxpol    = mxpol[~mxpol.mask]
    msigxpol = msigxpol[~msigxpol.mask]

    tmean = []
    xmean = []
    ymean = []
    tmin  = mtime[ignore]

    i = 0
    while mtime[i] < tmin :
        i += 1
    while tmin+step <= mtime[-1] :
        tdata,xdata,xsigma = [],[],[]
        while mtime[i] < tmin + step :
            tdata.append(mtime[i])
            xdata.append(mxpol[i])
            xsigma.append(msigxpol[i])
            i += 1
        tm = tmin + step/2.0
        xm = cm(tdata,xdata,xsigma,tmin,step)

        if not (math.isnan(tm) or math.isnan(xm)) :
            tmean.append(tm)
            xmean.append(xm)
        tmin += step
    return tmean,xmean


def plot_signalVSmean(ax,t,xdata,sigx,tm,xm,mask) :
    #simple plotting function
    step = tm[1] - tm[0]
    mt     = ma.masked_array(t,mask)
    mx     = ma.masked_array(xdata,mask)

#    xticks = tm
#    ax.set_xticks(xticks)
#    ax.errorbar(time,xdata,yerr=sigx,fmt='+',alpha=0.6,markersize=1)
#    ax.scatter(mt,mx,marker='o',s=1.,alpha=0.6)
#    ax.scatter(tm,xm,s=100,lw=2,marker='+',color='r',alpha=0.8)
    ax.plot(tm,xm,lw=1,alpha=0.8,label=str(step))


def interpolation_function(data_time,data_signal,data_sigma,interp_time_step=15):
#function to interpolate the data_signal and data_time with a step of inter_time_step

#argument of function#
#data_time        :: time of data set
#data_signal      :: signal of data set
#data_error       :: sigma of mesured signal
#interp_time_step :: time step of interpolation, 15 by default
######################

#return Value#
#interp_time   :: interpolated time with the good timestep
#interp_signal :: interpolated signal with the good timestep 
##############

#masking data to kill bad value off signal
  mask_value = get_mask(data_signal,data_sigma)

#processing interpolation
  interp_time,interp_signal=get_mean_signal(data_time, data_signal, data_sigma, mask_value, step=interp_time_step, ignore=200,method=2)

  return interp_time,interp_signal


#end function interpolation_function#
###################################################################

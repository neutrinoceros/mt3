import numpy as np
import pylab as pl
#from scipy.optimize import curve_fit as cf

def flin(x,a,b) :
    return a*x+b

def comp_mean(xdata,ydata,ysigma,tmin,h=7.0) :
    w = np.array(ysigma)**2
    for i in range(len(w)) :
        if w[i] == 0 :
            w[i] = 1.
    w = w**(-1)
    a,b = np.polyfit(xdata,ydata,1,w=w)
    ymean = a/2 * (h + 2*tmin) + b
    return ymean
    
filename = 'data/opa2015a.eops'
tab      = np.genfromtxt(filename,usecols=(0,4,5,9,10))

time,xpol,ypol,sigxpol,sigypol = tab[:,0],tab[:,1],tab[:,2],tab[:,3],tab[:,4]

tmean = []
xmean = []
ymean = []
step = 20.0 #jours
tmin = time[100]
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
    xm = comp_mean(tdata,xdata,xsigma,tmin,step)
    ym = comp_mean(tdata,ydata,ysigma,tmin,step)
    tmean.append(tm)
    xmean.append(xm)
    ymean.append(ym)
    tmin += step

#plotting
pl.ion()
fig,ax =  pl.subplots()

pl.show()
ax.scatter(time,xpol,marker='o',alpha=0.6)
#ax.errorbar(time,xpol,yerr=sigxpol,fmt='+',alpha=0.6)

xticks = [time[0] + i*step for i in range(len(tmean))]
ax.set_xticks(xticks)
ax.scatter(tmean,xmean,s=100,lw=2,marker='+',color='r',alpha=0.8)
pl.draw()
pl.ioff()

raw_input()

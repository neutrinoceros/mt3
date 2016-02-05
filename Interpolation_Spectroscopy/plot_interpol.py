#script in order to plot the different interpolation 
from interpol import *
#data loading
filename = '../data/opa2015a.eops'
tab      = np.genfromtxt(filename,usecols=(0,4,5,9,10))

time,xpol,ypol,sigxpol,sigypol = tab[:,0],tab[:,1],tab[:,2],tab[:,3],tab[:,4]

method   = raw_input("Use interpolation method 1 or 2 ? (I don't want to influence your choice but 2 is definitly better)     ")
while method not in ['1','2'] :
    method   = raw_input("""I said "1 or 2", that's not so hard a choice, I'll give you another try :      """)

method = int(method)


#plotting
pl.ion()
fig,axes = pl.subplots(nrows=2)
pl.show()

axes[0].set_title("method {}".format(method))

maskx      = get_mask(xpol,sigxpol)
masky      = get_mask(ypol,sigypol)

step = 15
#data processing

tmeanx,xmean = get_mean_signal(time,xpol,sigxpol,maskx,step=step,ignore=200,method=method)
tmeany,ymean = get_mean_signal(time,ypol,sigypol,masky,step=step,ignore=200,method=method)

plot_signalVSmean(axes[0],time,xpol,sigxpol,tmeanx,xmean,maskx)
plot_signalVSmean(axes[1],time,ypol,sigypol,tmeany,ymean,masky)


mtime     = ma.masked_array(time,maskx)
mxpol     = ma.masked_array(xpol,maskx)
axes[0].scatter(mtime,xpol,c='m',edgecolor='face',marker='o',s=1,alpha=0.8)
axes[0].errorbar(mtime,xpol,yerr=sigxpol,fmt='+',alpha=0.6,markersize=1)
pl.legend()
pl.draw()
pl.ioff()

raw_input()

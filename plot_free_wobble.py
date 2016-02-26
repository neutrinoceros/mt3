#-*-coding:utf-8-*-

#======================================================================
# script to 3D-plot the free_wobble
#======================================================================

import numpy as np
import pylab as pl
from astropy.time import Time
from mpl_toolkits.mplot3d import Axes3D

#data loading
#======================================================================

data = np.loadtxt('free_wobble.dat')
time, etax, etay = data[:,0],data[:,1],data[:,2]
def conv2years(t) :
    return 2000 - (51544.5-t)/365.25

ttime = np.array([conv2years(t) for t in time])
#ttime = Time(time,format='jd',scale='utc')
#print ttime.iso

#plot script
#======================================================================
pl.ion()
pl.show()


fig = pl.figure()
ax  = fig.gca(projection='3d')
ax.set_xlabel(r'$t$',size=20)
ax.set_ylabel(r'$x$',size=20)
ax.set_zlabel(r'$y$',size=20)
ax.set_xlim(1970,2020)
ax.set_ylim(-0.5,0.5)
ax.set_zlim(-0.5,0.5)
ax.set_xticks([1970+n*10 for n in range(6)])

ax.plot(ttime,etax,etay,c='royalblue',lw=1.5)
ax.plot(ttime,etay,zdir='z',c='crimson',zs=-0.5)
ax.plot(ttime,etax,zdir='y',c='darkcyan',zs=0.5)

pl.draw()
pl.savefig('pictures/free_wobble.pdf',transparent=True)
pl.ioff()
raw_input()

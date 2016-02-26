import numpy as np
import matplotlib.pyplot as plt 

correlation=np.loadtxt('fit_amplitude_correlation.dat')


plt.matshow(correlation,cmap=plt.cm.viridis)
axenumx=np.array([21,41.5,62,84,87])
axenumy=np.array([21,41.5,62,84,87])
axelegendx=('$A_{re}$','|','$A_{im}$','|','$cst$')
axelegendy=('$A_{re}$','-','$A_{im}$','-','$cst$')
plt.xticks(axenumx,axelegendx)
plt.yticks(axenumy,axelegendy)
plt.colorbar()

plt.savefig('pictures/fit_amplitude_correlation.pdf',tranparent=True)
# plt.show()

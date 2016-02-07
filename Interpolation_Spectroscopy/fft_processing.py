import numpy as np
import scipy as sc
import scipy.fftpack as fftp
import sys

#script in order to create fft (power and frequence)

def processing_fft(time,signal,sample_step):
######################
#argument of function#
######################
#time        : temporal data off the mesurment, the signal have to be sampled with a contant step
#signal      : value of the signal mesurment
#sample_step : sample_step of mesurment

#return variable#
######################
#frequences     : frequence of the fft
#spectral_power : spectral power (module of complex fft) of the mesuring series
######################

  #checking array dimension
  if (len(time)!=len(signal)):
    sys.exit("time and signal have note the same dimension")


  Nbr_of_point=len(time)

  #spectral power (normalized) computation
  Cplx_ft=fftp.fft(signal)
  spectral_power=np.abs(Cplx_ft[:Nbr_of_point/2])*2/Nbr_of_point

  #frequences computation
  frequences=np.linspace(0.,1./(2.*sample_step),Nbr_of_point/2)

  return frequences,spectral_power

#end of processing_fft#
############################################

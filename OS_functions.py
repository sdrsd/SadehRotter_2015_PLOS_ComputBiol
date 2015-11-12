
######################################################################################
# OS_functions -- Functions needed for analysis and plotting of the results
#
# Reference: Sadeh and Rotter 2015.
# "Orientation selectivity in inhibition-dominated networks of spiking neurons:
# effect of single neuron properties and network dynamics" PLOS Computational Biology.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
######################################################################################

import scipy.optimize
import numpy as np

###
# Fourier components
def _fft_(x):
    n = len(x)
    z = np.fft.fft(x)
    f0 = abs(z[0])/n      # mean     
    f1 = abs(z[1])/(n/2) # 1st order
    return f0, f1

###
# OSI: 1- CV; CV: sum(r exp(j*2*t))/sum(r): t=0...pi (rad)
def _osv_(r, t):
    z = np.sum(r * np.exp(1j *2.*t))/np.sum(r)
    osi = abs(z)  #osi
    po = np.angle(z) %(2*np.pi) #PO
    return osi, po

###
# fit von-Mises
def vonMises(xdata,ydata):
    # define a von-Mises fitting function: B + A exp(k cos(2(theta - phi)) -1)
    # p[0] = B
    # p[1] = A
    # p[2] = k (1/D)
    # p[3] = phi
    fitfunc = lambda p, x: p[0] + p[1] * \
                           np.exp( p[2]*( np.cos(2*(x - p[3]))-1.) )
    errfunc = lambda p, x, y: fitfunc(p,x)-y
    
    # initial guess
    p0 = [1., 1., 1., 1.]
    # fit
    p1, success = scipy.optimize.leastsq(errfunc, p0[:],args=(xdata,ydata))
    
    tw = (90./np.pi) * np.arccos(abs(1 + np.log((1.+np.exp(-2.*p1[2]))/2)/p1[2]) )
    
    return p1, fitfunc(p1,xdata), tw, success

######################################################################################

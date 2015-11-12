
######################################################################################
# OS_params -- Default parameters for simulations in the following:
#
# Reference: Sadeh and Rotter 2015.
# "Orientation selectivity in inhibition-dominated networks of spiking neurons:
# effect of single neuron properties and network dynamics" PLOS Computational Biology.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
######################################################################################

import numpy as np
import pylab as pl
import time, os
import pickle as cPickle
from imp import reload

code_path = os.getcwd()
if 'Results' not in os.listdir(code_path):
    os.mkdir(code_path + '/Results/')
res_path = code_path+ '/Results/'

# define your nest path here, if it was needed
#nest_path = code_path
nest_path = '/Users/sadra/NEST/nest/ins/lib/python3.4/site-packages/'
os.chdir(nest_path); import nest; os.chdir(code_path)

## fixing the seed
np.random.seed(1234)
  
################################################################################

### neurons and connectivity
N = 5000  # number of neurons in the local network
exc_inh = .8
NE = int(exc_inh *N)
NI = N - NE

### recurrent coupling
g = 4.
eps_exc = .2
eps_inh = .5
j_rec = .1 # local EPSP
J_rec = [j_rec, j_rec, -g*j_rec, -g*j_rec]

### non-local input
n_ext = 5000. # number of neurons in the non-local network
r_ext = n_ext * 1. # spont. firing rate of the non-local network
j_ext = .2   # non-local EPSP

### ffw input
n_ffw = 50. # number of LGN neurons contacting one cortical neuron
r_ffw = 20. # avg. firing rate of LGN neurons
sb = n_ffw * r_ffw  # ffw overall rate
j_ffw = 1.  # ffw EPSP

### modulations
n_smpl = 12 # sample neurons for recording the membrane potential
m = 0.2    # input modulation

### Initial random PO of neurons
po_init = np.random.uniform(0., np.pi, N) 
po_init[0:n_smpl] = np.arange(0, 180, 180/12.)*np.pi/180
po_init[NE:NE+n_smpl] = np.arange(0, 180, 180/12.)*np.pi/180


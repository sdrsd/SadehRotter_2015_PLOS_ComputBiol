
######################################################################################
# OS_run -- Reads the functions from OS_source and runs the network simulations
#
# Reference: Sadeh and Rotter 2015.
# "Orientation selectivity in inhibition-dominated networks of spiking neurons:
# effect of single neuron properties and network dynamics" PLOS Computational Biology.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
######################################################################################

from imp import reload
import OS_params; reload(OS_params); from OS_params import *

import OS_source; reload(OS_source)

#####################################
########## parameters of Figure 1
########## change it for other figures
g = 4
neuron_type = 'pif' # 'lif'
delay_type = 'random' # 'fixed'
#####################################

### instantiation of the class
my_net = OS_source.Network_Stim(N = N, exc_inh = exc_inh, g = g, J_ext = j_ext, J_ffw = j_ffw, J_rec = J_rec, \
                               r_base = sb, po_init = po_init, t_trans=150., r_ext=r_ext, m=m, eps_exc=eps_exc, eps_inh=eps_inh)

### delays (fixed or random?)
my_net.delay_fixed = 1.
my_net.delay_ext = 1. 
Nd = eps_exc*NE + eps_inh*NI

if delay_type == 'fixed':
   my_net.delay_dist = 3.*np.ones((N,Nd))
elif delay_type == 'random':
   my_net.delay_dist = np.random.uniform(0.1, 3., (N, Nd))

stim_no = 8

### lif or pif
if neuron_type == 'lif':    
   ## lif
   my_net.tauMem = 20.
elif neuron_type == 'pif':
   ## pif
   my_net.tauMem = 1e12

# results are saved in this folder
sim_folder = 'N-'+str(my_net.N)+'_'+neuron_type + '_delayType-'+delay_type + '_g-'+str(my_net.g)

### make the connectivity (1) or read it (0)?
con_gen = 1 #0

if con_gen == 1:
	### generate and write
	con_exc, con_inh = my_net.random_connectivity(eps_exc = eps_exc, eps_inh = eps_inh)

	con = {}
	con['exc'], con['inh'] = con_exc, con_inh

	fl = open('con', 'wb')
	cPickle.dump(con, fl, 2)
	fl.close()
else:
	### read      
	fl = open('con', 'rb')
	con = cPickle.load(fl)
	fl.close()

con_exc, con_inh = con['exc'], con['inh']

### simulation folder
if sim_folder in os.listdir(res_path):
   os.chdir(res_path+sim_folder)
else:
    os.mkdir(res_path+sim_folder) 
    os.chdir(res_path+sim_folder)

### simulation parameters
my_net.contrast = [1., 2., 3.]
my_net.trial_no = len(my_net.contrast)

my_net.stim_range = np.arange(0., np.pi, np.pi/stim_no)

# simulation time for each orientation (in ms)
my_net.simtime = 3.*1000. + my_net.t_trans 

### spontaneous activity
print('##############################')
print('### spontaneous activity')
print('##############################') 
my_net.stim_rndm(stim_deg=0., con_exc=con_exc, con_inh=con_inh, stim_id = 'spont', \
                     mem_pot_rec = 1, n_smpl=n_smpl, contrast = [0.,0.,0.])

### stimulation (evoked)
my_net.input_modulation = m
for stim in enumerate(my_net.stim_range):
    print('##############################')
    print('### orientation: ', \
    str(stim[0]+1),'/',str(stim_no)+' ('+str(stim[1]*180./np.pi)+' deg)')
    print('##############################')
    stim_id = 'st'+str(stim[0])
    my_net.stim_rndm(stim_deg=stim[1], con_exc=con_exc, con_inh=con_inh, stim_id = stim_id, \
                     mem_pot_rec = 1, n_smpl=n_smpl, contrast = my_net.contrast)


my_net.save_info()


os.chdir(code_path)

################################################################################



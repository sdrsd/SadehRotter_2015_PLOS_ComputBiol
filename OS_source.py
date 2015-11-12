
######################################################################################
# OS_source -- The source functions to connect and run the network simulations in NEST
#
# Reference: Sadeh and Rotter 2015.
# "Orientation selectivity in inhibition-dominated networks of spiking neurons:
# effect of single neuron properties and network dynamics" PLOS Computational Biology.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
######################################################################################

from imp import reload
import OS_params; reload(OS_params); from OS_params import *

################################################################################
## Main Class
################################################################################
class Network_Stim (object):
      def __init__(self, N=5000, exc_inh=.8, g=8., J_rec=[.1, .1, .1, .1], J_ext=.1, J_ffw=1., \
                   r_base=15000., r_ext=0., m=.1, po_init=[], t_trans=300., eps_exc=.2, eps_inh=.5):
          
          ##############################################################
          ## Params
          ##############################################################
          self.simtime = 2.*1500.0 # Simulation time in ms 
          self.t_trans = t_trans
          
          self.g       = g
          self.eps_exc = eps_exc   # connection probability
          self.eps_inh = eps_inh
          
          self.sig_exc = eps_exc
          self.sig_inh = eps_inh
          self.typ = []
          
          self.N       = N
          self.exc_inh = exc_inh
          self.NE      = int(self.exc_inh*self.N)
          self.NI      = self.N - self.NE          
          
          self.dt      = 0.1    # the resolution in ms
          self.delay_fixed   = 1.5    # synaptic delay in ms
          self.delay_dist = .1
          self.delay_ext = .1

          self.tauSyn  = 0.5
          self.tauMem  = 20.0
          self.theta   = 20.0
          
          self.J_rec   = J_rec
          self.J_ee    = self.J_rec[0]
          self.J_ei    = self.J_rec[1]
          self.J_ie    = self.J_rec[2]
          self.J_ii    = self.J_rec[3]
          
          self.J_ext   = J_ext          
          self.J_ffw   = J_ffw

          self.po_init = po_init              
          
          self.input_modulation = m
          self.r_base = r_base
          self.r_ext = r_ext
          self.stim_deg = 0.  
          self.stim_range = [0.]  
          self.trial_no = 1      
          self.n_smpl = 12
          self.contrast = []

################################################################################          
## Functions
################################################################################
      
      # info
      def save_info(self):
          info ={}
          
          info['tauMem'] = self.tauMem
          info['theta']  =self.theta
          info['dt'] = self.dt
          info['delay'] = self.delay_fixed

          info['N'] = self.N
          info['g'] = self.g
          info['J'] = self.J_rec
          info['b'] = self.r_base
          info['e'] = self.r_ext
          
          info['eps_exc'] = self.eps_exc
          info['eps_inh'] = self.eps_inh
          
          info['exc_inh'] = self.exc_inh
          info['m'] = self.input_modulation         
          info['J_ext'] = self.J_ext
          info['J_ffw'] = self.J_ffw
           
          info['simtime']=self.simtime
          info['t_trans']=self.t_trans
          
          info['stim_range'] = self.stim_range
          info['trial_no'] = self.trial_no
          info['po_init'] = self.po_init
          info['n_smpl'] = self.n_smpl
          info['contrast'] = self.contrast

          fl = open('info', 'wb')
          cPickle.dump(info, fl, 2)
          fl.close()
      
      # random connectivity
      def random_connectivity(self, eps_exc=0.1, eps_inh=0.1):
          
          print("Connectivitiy ...")
          
          self.eps_exc = eps_exc
          self.eps_inh = eps_inh
          self.CE      = int(self.eps_exc*self.NE)   # exc synapses per neuron
          self.CI      = int(self.eps_inh*self.NI)   # inh synapses per neuron
          
          np.random.seed(1234)
          
          ### fixed in-degree
          Cin_exc, Cin_inh = [], []
          for n in range(1, self.N+1):
              if n%1000 == 0:
                 print(n)
              # exc
              if n <= self.NE:
                 zex = list(range(1, n)) + list(range(n+1,self.NE+1))
              else:
                 zex = list(range(1, self.NE+1))
              np.random.shuffle(zex)
              cin_exc = zex[0:self.CE]
              
              # inh
              if n <= self.NE:
                 zin = list(range(self.NE+1, self.N+1))
              else:
                 zin = list(range(self.NE+1, n)) + list(range(n+1,self.N+1))
              np.random.shuffle(zin)
              cin_inh = zin[0:self.CI]
              
              Cin_exc.append(cin_exc)
              Cin_inh.append(cin_inh)
          
          return Cin_exc, Cin_inh
          
      # stimulation
      def stim_rndm(self, stim_deg, con_exc=[], con_inh=[], stim_id=0, mem_pot_rec=1, n_smpl=12, contrast=[1.]):
          self.n_smpl = n_smpl
          self.contrast = contrast

          ### Start
          ti = time.time()
          
          self.stim_deg = stim_deg
          
          ### Nest start
          nest.ResetKernel()
          nest.SetStatus([0],[{'rng_seeds':[1312738]}])
          nest.SetKernelStatus({"resolution": self.dt, 
                                "print_time": True, 
                                "overwrite_files": True, 
                                'local_num_threads': 1})
          nest.sr("/RandomConvergentConnect << /allow_multapses false >> SetOptions")
          
          print("Network ...")
                    
          ### current-based delta
          neuron_params_delta = {"C_m"       : 1.0,
                          "tau_m"     : self.tauMem,
                          "t_ref"     : 2.0,
                          "E_L"       : 0.0,
                          "V_min"     : - np.inf,
                          "V_m"        : 0.,
                          "V_reset"   : 0.,
                          "V_th"      : self.theta,
                          "I_e"       : 0.}                         
          
          nest.SetDefaults("iaf_psc_delta", neuron_params_delta)
          nodes_ex = nest.Create("iaf_psc_delta", self.NE)
          nodes_in = nest.Create("iaf_psc_delta", self.NI)
                    
          nodes_all = nodes_ex + nodes_in

          ##
                          
          self.input_rate = self.r_base *(1.+self.input_modulation* np.cos(2*(self.stim_deg - self.po_init)) )
          
          noise = nest.Create("poisson_generator", self.N) 
          
          # devices: spike detector
          sp_all_trans = nest.Create("spike_detector", self.trial_no)  
          sp_all = nest.Create("spike_detector", self.trial_no)      
          for trial in range(self.trial_no):
              ### Spike Detector
              # transient     
              nest.SetStatus([sp_all_trans[trial]], 
                          {"label":"spikes-all-trans-"+stim_id+'-tr'+str(trial), 
                           "withgid":True, 
                           "withtime":True,
                           "to_file":True,
                           "to_memory":False,
                           "start": self.simtime*trial + 0.,
                           "stop": self.simtime*trial + self.t_trans })
              nest.ConvergentConnect(nodes_all, [sp_all_trans[trial]])
              # stationary
              nest.SetStatus([sp_all[trial]], 
                        {"label":"spikes-all-"+stim_id+'-tr'+str(trial), 
                         "withgid":True, 
                         "withtime":True,
                         "to_file":True,
                         "to_memory":False,
                         "start": self.simtime*trial + self.t_trans,
                         "stop": self.simtime*trial + self.simtime }) 
              nest.ConvergentConnect(nodes_all, [sp_all[trial]])
          
          # devices: volt meter     
          ### Voltmeter (one exc & one inh)                    
          if mem_pot_rec == 1:
             # stationary
             vm_ex_stat = nest.Create("voltmeter")
             nest.SetStatus(vm_ex_stat, {"withtime":True, 
                         "label":"vm-exc-"+stim_id, #+'-tr'+str(trial),
                         "to_file":True,
                         "to_memory":False,
                         "start": 0., #self.simtime*trial + self.t_trans,
                         "stop": np.inf}) #self.simtime*trial + self.simtime})          
             nest.DivergentConnect(vm_ex_stat, nodes_ex[0:self.n_smpl])
             
             vm_in_stat = nest.Create("voltmeter")
             nest.SetStatus(vm_in_stat, {"withtime":True, 
                         "label":"vm-inh-"+stim_id, #+'-tr'+str(trial),
                         "to_file":True,
                         "to_memory":False,
                         "start": 0., #self.simtime*trial + self.t_trans,
                         "stop": np.inf}) #self.simtime*trial + self.simtime})          
             nest.DivergentConnect(vm_in_stat, nodes_in[0:self.n_smpl])
                 
                   
          ### Connecting
          print("Connections ...")       
          
          for nn in enumerate(nodes_all):
              nest.ConvergentConnect([noise[nn[0]]], [nn[1]], self.J_ffw, self.delay_fixed)
          
          ext_inp = nest.Create("poisson_generator")  
          nest.SetStatus(ext_inp, {"rate":self.r_ext})
          nest.DivergentConnect(ext_inp, nodes_all, self.J_ext, self.delay_ext)
          
          ### recurrent connectivity
          zzz = self.eps_exc*self.NE
          for ne in range(self.NE):
              # exc --> exc
              nest.ConvergentConnect(con_exc[ne], [nodes_all[ne]], self.J_ee, self.delay_dist[ne,0:zzz].tolist())
              # inh --> exc
              nest.ConvergentConnect(con_inh[ne], [nodes_all[ne]], self.J_ie, self.delay_dist[ne, zzz:].tolist())
                        
          for ni in range(self.NE, self.N):
              # exc --> inh
              nest.ConvergentConnect(con_exc[ni], [nodes_all[ni]], self.J_ei, self.delay_dist[ni,0:zzz].tolist())              
              # inh --> inh
              nest.ConvergentConnect(con_inh[ni], [nodes_all[ni]], self.J_ii, self.delay_dist[ni,zzz:].tolist())              
          
          print("Simulation ...")
          for it in range(self.trial_no):
              print('####################')
              print('### Contrast ', str(it+1))
              print('####################')
              ti = time.time()
              for nn in enumerate(nodes_all):
                  nest.SetStatus([noise[nn[0]]], {'rate':self.contrast[it]*self.input_rate[nn[0]]})             
                        
              nest.Simulate(self.simtime)     
                        
              ts = time.time()
              sim_time = ts - ti
              
              print('########################################')
              print("Simulation time   : %.2f s" % sim_time)
              print('########################################')

################################################################################


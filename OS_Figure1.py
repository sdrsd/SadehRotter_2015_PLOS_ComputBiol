
######################################################################################
# OS_Figure1 -- Reads and analyzes the results generated from OS_run.py
# and plots Figure 1 in the following:
#
# Reference: Sadeh and Rotter 2015.
# "Orientation selectivity in inhibition-dominated networks of spiking neurons:
# effect of single neuron properties and network dynamics" PLOS Computational Biology.
#
# Author: Sadra Sadeh <s.sadeh@ucl.ac.uk> // Created: 2014-2015
######################################################################################

from imp import reload
import OS_params; reload(OS_params); from OS_params import *

import OS_functions; reload(OS_functions); from OS_functions import *

from mpl_toolkits.axes_grid.inset_locator import inset_axes

cwd = code_path

################################################################################
################################################################################

## for the first time you run the file, you need to make these data (xxx_make = 1)
## for the subsequent runs, you can turn the flag off, the data is now loaded
spike_data_make = 1#0
isi_cv_make = 1#0

#################################

## simulation folder for the results
sim_folder = 'N-5000_pif_delayType-random_g-4'

## figure name
fig_name = 'Figure1'

################################################################################
################################################################################

ti = time.time()   

#######################################
### reading simulation params & results

os.chdir(res_path+sim_folder)

fl = open('info', 'rb')
infor = cPickle.load(fl)
fl.close()

po_init = infor['po_init']
inputs = [infor['b']]

print(infor)

# params
N = infor['N']
NE = int(.8*N)
NI = N - NE

stim_range = infor['stim_range'] 
stim_range_deg = stim_range*180./np.pi 
simtime = infor['simtime']
t_trans = infor['t_trans']
tr_time = simtime - t_trans

CE = infor['eps_exc']
CI = infor['eps_inh']
f = infor['exc_inh']
g = infor['g']

tau_m = infor['tauMem']/1000.
V_th = infor['theta']
V_r = 0.
t_ref = 2./1000

Js = infor['J_ffw']
J_ext = infor['J_ext']
sb = infor['b']
se = infor['e']
J = infor['J'][0]
m = infor['m'] 
sm = m*sb
trial_no = infor['trial_no']
tr_no = infor['trial_no']
contrast = infor['contrast']

n_smpl = infor['n_smpl']

fl = open('results', 'rb')
results = cPickle.load(fl)
fl.close()

# results
tc_trans = results['tc_trans']
tc = results['tc']

tc_mean = np.mean(tc, 1)
tc_std = np.std(tc, 1)

tc_f0 = results['tc_f0']
tc_f1 = results['tc_f1']

vm_tc = results['vm_tc']

VM_fit = results['VM_fit']
err_fit = results['err_fit']
scs_fit = results['scs_fit']
TW_out = results['TW_out']
PO_out = results['PO_out']
OSI_out = results['OSI_out']

################################################################################
os.chdir(code_path)

################################################################################
### spike data all
################################################################################

## spont
os.chdir(res_path+sim_folder)

spd_all_sp = {}
for st in range(len(stim_range)):
    #spd_all[st] = []
    for tr in range(tr_no): 
        spd_all_sp[tr] = []       
        spflno = 2*N+ 1+tr + tr_no
        spike_file = 'spikes-all-spont-tr'+str(tr)+'-'+str(spflno)+'-0.gdf'
        spd_tr = pl.loadtxt(spike_file)
        #spd += spd_tr.tolist()
        spd_all_sp[tr] = spd_tr

os.chdir(cwd)	

if spike_data_make == 1:
    os.chdir(res_path+sim_folder)
    
    spd_all = {}
    for st in range(len(stim_range)):
        spd_all[st] = []
        for tr in range(tr_no):        
            spflno = 2*N + 1+tr + tr_no
            spike_file = 'spikes-all-st'+str(st)+'-tr'+str(tr)+'-'+str(spflno)+'-0.gdf'
            spd_tr = pl.loadtxt(spike_file)
            #spd += spd_tr.tolist()
            spd_all[st].append(spd_tr)
    
    os.chdir(cwd)	
    
    fl = open('spd_all', 'wb')
    cPickle.dump(spd_all, fl)
    fl.close()
else:
    os.chdir(cwd)
    fl = open('spd_all', 'rb')
    spd_all = cPickle.load(fl)
    fl.close()

os.chdir(cwd)

################################################################################
### ISI CV
################################################################################

if isi_cv_make == 1:
    print('### ISI CV ###')
    ### ISI CV
    #ISI, rr = {}, {}
    CV = {}
    for tr in range(trial_no):
        print('tr: '+str(tr))
        CV[tr] = []
        for i in range(1, N+1):
            spt = spd_all[4][tr][:,1][np.where(spd_all[4][tr][:,0] == i)]
            isi = np.diff(spt)
            cv = np.std(isi)/np.mean(isi)
            if not np.isnan(cv) and not np.isinf(cv) and len(spt) > 10:
               CV[tr].append(cv)
        CV[tr] = np.array(CV[tr])

    isi_cv = CV
    
    fl = open('isi_cv', 'wb')
    cPickle.dump(isi_cv, fl)
    fl.close()
else:
    fl = open('isi_cv', 'rb')
    isi_cv = cPickle.load(fl)
    fl.close()

################################################################################
################ raster plots for population
################################################################################

al = [.25, .5, 1.]
fg_lb_sz = 20

pl.figure(figsize= (15,9))

########## ########## ########## ########## ########## 
########## raster plot for the middle contrast (spont)
########## ########## ########## ########## ########## 

tr = 1

zz_sp = spd_all_sp[tr]
    
excid_sp = np.where(zz_sp[:,0] <= NE)[0]
inhid_sp = np.where(zz_sp[:,0] > NE)[0]

ax = pl.subplot(3, 3, 1)
pl.title(r'Spont. spiking activity, contrast C = 2', size=15) 

ax.text(-.15, 1., 'A1', size=fg_lb_sz, transform = ax.transAxes)

ax.plot(zz_sp[excid_sp, 1], zz_sp[excid_sp, 0], 'r.', ms=2)
ax.plot(zz_sp[inhid_sp, 1], zz_sp[inhid_sp, 0], 'b.', ms=2)
#
t1 = tr*(tr_time+t_trans) + 500
dt = 60
t2 = t1 + dt 

xtk = np.array([t1, t1+20, t1+40, t2])
ax.set_xticks(xtk)
ax.set_xticklabels(xtk - simtime)

ax.set_yticks([1, 1000, 2000, 3000, 4000, 5000])
ax.set_yticklabels([1, '', '', '', '4k', '5k'])

ax.set_xlim(t1-5, t2+5)
ax.set_ylim(0-100, N + 100)

#ax.set_xlabel('time (ms)')
ax.set_ylabel('Neuron #')

########## ########## ########## ########## ########## 
ax = pl.subplot(6, 3, 2*3+1)

binw = 5.
bins = simtime/binw

ax.hist(zz_sp[excid_sp, 1], bins = bins, weights=1000./binw*np.ones(len(zz_sp[excid_sp, 1]))/NE, color='r', alpha=.75)
ax.hist(zz_sp[inhid_sp, 1], bins = bins, histtype='step', weights=1000./binw*np.ones(len(zz_sp[inhid_sp, 1]))/NI, color='b', lw=2, alpha=1)

ax.set_xticks(xtk)
ax.set_xticklabels(xtk - simtime)
ax.set_xlim(t1-5, t2+5)

ax.set_ylim(0-1, 40+1)

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Population rate \n (spikes/sec)', size=10)

ax.text(.7, .8, 'binw: '+str(int(binw))+' ms', transform = ax.transAxes)

########## ########## ########## ########## ########## 
########## raster plot for the middle contrast (evoked)
########## ########## ########## ########## ########## 

tr = 1

zz = spd_all[4][tr]
    
excid = np.where(zz[:,0] <= NE)[0]
inhid = np.where(zz[:,0] > NE)[0]

ax = pl.subplot(3, 3, 2)
pl.title(r'Evoked spiking activity, stim. orient. $\theta: 0^\circ$, C = 2', size=12.5) 

ax.text(-.15, 1., 'A2', size=fg_lb_sz, transform = ax.transAxes)

ax.plot(zz[excid, 1], zz[excid, 0], 'r.', ms=2)
ax.plot(zz[inhid, 1], zz[inhid, 0], 'b.', ms=2)

xtk = np.array([t1, t1+20, t1+40, t2])
ax.set_xticks(xtk)
ax.set_xticklabels(xtk - simtime)

ax.set_yticks([1, 1000, 2000, 3000, 4000, 5000])
ax.set_yticklabels([1, '', '', '', '4k', '5k'])

ax.set_xlim(t1-5, t2+5)
ax.set_ylim(0-100, N + 100)

#ax.set_xlabel('time (ms)')
ax.set_ylabel('Neuron #')

########## ########## ########## ########## ########## 
ax = pl.subplot(6, 3, 2*3+2)

binw = 5.
bins = simtime/binw

ax.hist(zz[excid, 1], bins = bins, weights=1000./binw*np.ones(len(zz[excid, 1]))/NE, color='r', alpha=.75)
ax.hist(zz[inhid, 1], bins = bins, histtype='step', weights=1000./binw*np.ones(len(zz[inhid, 1]))/NI, color='b', lw=2, alpha=1)

ax.set_xticks(xtk)
ax.set_xticklabels(xtk - simtime)

#ax.set_xlim(0-10, tr_time+10)
ax.set_xlim(t1-5, t2+5)

ax.set_ylim(0-1, 40+1)

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Population rate \n (spikes/sec)', size=10)

ax.text(.7, .8, 'binw: '+str(int(binw))+' ms', transform = ax.transAxes)

########## ########## ########## ########## ########## 
########## sorted raster plots for all contrasts
########## ########## ########## ########## ########## 

for tr in range(trial_no):
    zz = spd_all[4][tr]
    
    excid = np.where(zz[:,0] <= NE)[0]
    inhid = np.where(zz[:,0] > NE)[0]

    #
    t1 = tr*(tr_time+t_trans) + 500
    dt = 60
    t2 = t1 + dt 

    ###
    ax = pl.subplot(6, 3, 3+3*tr)
    #pl.title(r'$C = $'+str(contrast[tr])) 
    ax.text(.05, 1.05, 'C = '+str(int(contrast[tr])), transform = ax.transAxes, fontsize=15)#, fontstyle='bold')
    if tr == 0: ax.text(-.1, 1., 'B', size=fg_lb_sz, transform = ax.transAxes)

    ax.plot(zz[excid, 1], po_init[(zz[excid, 0]-1).astype('int')], 'r.', ms=2.)
    ax.plot(zz[inhid, 1], po_init[(zz[inhid, 0]-1).astype('int')], 'b.', ms=2.)
    
    ax.set_xticks(xtk)
    ax.set_xticklabels([])#xtk - tr*(simtime - t_trans))
    ax.set_yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.set_yticklabels([])

    ax.set_xlim(t1-5, t2+5)
    ax.set_ylim(0-.1, np.pi+.1)    

    if tr == 1: 
       ax.set_yticklabels([0, 45, 90, 135, 180])
       ax.set_ylabel('Input PO (deg)')

xtk = np.array([t1, t1+20, t1+40, t2])

ax.set_xticks(xtk)
ax.set_xticklabels(xtk - tr*(simtime))

ax.set_xlabel('Time (ms)')

########## ########## ########## ########## ########## 
########## dist. of firing rates
########## ########## ########## ########## ########## 

ax = pl.subplot(2, 4, 5)
#adjust_spines(ax,['left', 'bottom'], outward=0, s=0)
ax.text(-.2, .975, 'C', size=fg_lb_sz, transform = ax.transAxes)

for tr in range(trial_no):
    ax.hist(tc[4, tr, :], 50, histtype = 'step', color='k', lw=2, alpha=al[tr], label=str(int(contrast[tr])))

ax.set_yscale('log')

pl.legend(title = 'C', loc=2, frameon=False, prop={'size':12})

ax.set_ylim(0, 10000)
ax.set_xticklabels([])
ax.set_ylabel('#')

ax.set_xticks([0, 20, 40, 60, 80])
ax.set_xticklabels([0, 20, 40, 60, 80])

ax.set_xlabel('Firing rate (spikes/sec)')
#ia.set_ylabel('#')

########## dist. of ISI CV
########## ########## ########## ########## ########## 
### isi cv in the inset
########## ########## ########## ########## ########## 
#ax = pl.subplot(2, 4, 6)
ia = inset_axes(ax, width="30%", height="30%", loc=1)
#adjust_spines(ia,['left', 'bottom'], outward=0, s=0)

for tr in range(trial_no):
    ia.hist(isi_cv[tr], 50, histtype = 'step', color='k', lw=3, alpha=al[tr], label=str(int(contrast[tr])))

#pl.legend(title = 'C', )

ia.set_yticks([0, 100, 200, 300, 400, 500, 600])
ia.set_yticklabels([0, '', '', '', '', 500, ''])
ia.set_ylim(0-5, 600)

ia.set_xticks([0, .5, 1, 1.5, 2])
ia.set_xticklabels([0, '', 1, '', 2])

ia.set_xlim(0, 2.)

ia.set_xlabel('CV[ISI]')
ia.set_ylabel('#')


########## ########## ########## ########## ########## 
########## population tuning curves
########## ########## ########## ########## ########## 

#al = [.25, .5, 1.]
ax = pl.subplot(2, 4, 6)
#adjust_spines(ax,['left', 'bottom'], outward=0, s=0)
#tr = 1
ax.text(-.15, .975, 'D', size=fg_lb_sz, transform = ax.transAxes)

pl.title('Network output tuning curve')

ax.text(.1, .9, 'C = 2', size=15, transform = ax.transAxes)

ax.plot(po_init[0:NE], tc[4, tr, 0:NE], 'ro', ms=3.5, alpha=1, label='Exc.')
ax.plot(po_init[NE:], tc[4, tr, NE:], 'bo', ms=3.5, alpha=1, label='Inh.')

inp_mean = np.mean(tc[4,tr,:])
po_rng = np.arange(0, np.pi, .1)
ax.plot(po_rng, inp_mean*(1 + m*np.cos(2*(np.pi/2-po_rng))), 'g-', lw=2, label='Inp.')

pl.legend(loc=1, frameon=False, numpoints=1, markerscale=1.5, prop={'size':12.5}) 

ax.set_xlabel('Input PO (deg)')
ax.set_ylabel('Firing rate (spikes/sec)')

ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels([0, 45, 90, 135, 180])

ax.set_xlim(0-.1, np.pi+.1)
ax.set_ylim(0, 70+1)

########## ########## ########## ########## ########## 
########## mean population tuning curves
########## ########## ########## ########## ########## 

ax = pl.subplot(2, 4, 7)
#adjust_spines(ax,['left', 'bottom'], outward=0, s=0)
ax.text(-.15, .975, 'E', size=fg_lb_sz, transform = ax.transAxes)

pl.title('Avg. output tuning curve')

dxx = np.pi/100.
xx = np.arange(0, np.pi, dxx)

for tr in range(trial_no):
    zz = tc[4, tr, :]     
    zz_mean = np.zeros(len(xx))
    zz_std = np.zeros(len(xx))
    for ii, xi in enumerate(xx):
        ids = np.where( (po_init > xi) * (po_init < xi+dxx) == True)
        zz_mean[ii] = np.mean(zz[ids])
        zz_std[ii] = np.std(zz[ids])

    p, ft, tw, sc = OS_functions.vonMises(xx, zz_mean)
    print(tw)

    ax.plot(xx, zz_mean, 'k-', lw=2, alpha = al[tr], label=str(np.round(tw,1))+r'$^\circ$')
    ax.fill_between(xx, zz_mean - zz_std, zz_mean + zz_std, color='k', alpha=.1)

    inp_mean = np.mean(zz_mean)
    ax.plot(xx, inp_mean*(1 + m*np.cos(2*(np.pi/2-xx))), 'g-', lw=2, alpha = al[tr])#, label='Inp.')

pl.legend(frameon=False, title='output TW', prop={'size':12.5})

ax.set_xlabel('Input PO (deg)')

ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels([0, 45, 90, 135, 180])

ax.set_xlim(0-.1, np.pi+.1)
ax.set_ylim(0, 60+1)


########## ########## ########## ########## ########## 
########## linear prediction
########## ########## ########## ########## ########## 

os.chdir(code_path)

inv_make = 1

fl = open('con', 'rb')
con_loc = cPickle.load(fl)
fl.close()

con_exc, con_inh = np.array(con_loc['exc']), np.array(con_loc['inh'])

##############

print('### building weight matrix W ###')

Je = J
Ji = -g*J

W = np.zeros((N,N))
for ne in range(NE):    
    W[ne][con_exc[ne]-1] = Je 
    W[ne][con_inh[ne]-1] = Ji 
for ni in range(NE,N):    
    W[ni][con_exc[ni]-1] = Je 
    W[ni][con_inh[ni]-1] = Ji 

##############

if inv_make == 1:     
   print('### computing the inverse A = (1- W/Vth)^-1 ###')
   
   A = np.linalg.inv(np.eye(N) - W/V_th)
   
   fl = open('A', 'wb')
   cPickle.dump(A, fl)
   fl.close()

elif inv_read == 1:
   print('### reading the inverse A = (1- W/Vth)^-1 ###')
   
   fl = open('A', 'rb')
   A = cPickle.load(fl)
   fl.close()

################################################################################
print('### computing the linear prediction r = A s ###')
################################################################################

def _rect_(zz): return zz*(zz>0)

####
tc_out_L = []
for i in range(len(stim_range)):
    
    tc_out = []
    #tc_out_mod = []
    
    for ic, ct in enumerate(contrast):
        inp_b = ct * sb * Js/V_th
        inp_e = se * J_ext/V_th
        out_b = np.sum(A[0]* (inp_b+inp_e) )
        
        inp_m = ct * sm*np.cos(2*(stim_range[i]-po_init)) * Js/V_th      
        
        inp_tot = inp_b + inp_e + inp_m       
 
        r = np.array( (np.matrix(A) * np.matrix(inp_tot).T).T)[0]
        
        tc_out.append(r)    
    tc_out_L.append(tc_out)
tc_out_L = np.array(tc_out_L)

tc_out_LR = _rect_(tc_out_L)

######### plot it

ax = pl.subplot(2, 4, 8)
ax.text(-.15, .975, 'F', size=fg_lb_sz, transform = ax.transAxes)

pl.title('Predicting network responses')

tc_out_LR_rp = tc_out_LR/(1. + tc_out_LR*t_ref)

ax.plot(tc[4, 1, :], tc_out_LR_rp[4, 1, :], 'k.')#, ms=1)

ax.plot([0, 60], [0, 60], 'k--')

ax.set_xlabel('Simulated firing rate (spikes/sec)')
ax.set_ylabel('Linear prediction (spikes/sec)')

ia = inset_axes(ax, width="30%", height="30%", loc=2)
ia.yaxis.set_label_position('right')
ia.yaxis.tick_right()

for i, ct in enumerate(contrast):
    ia.hist(tc[4, i, :] - tc_out_LR_rp[4, i, :], bins=100, histtype='step', alpha=al[i], color='k')

ia.set_yscale('log')
ia.set_xticks([-10, -5, 0, 5, 10])
ia.set_xticklabels(['', -5, 0, 5, ''])
ia.set_xlim(-10, 10)
ia.set_ylabel('#')
ia.set_xlabel('Rate Diff')

#######################

pl.subplots_adjust(left=.05, right=.975, bottom=.075, top=.95, wspace=.25, hspace=.45)

pl.savefig(fig_name+'.png')


pl.show()

tf = time.time()   

print('time: ', np.round((tf - ti)/60., 1), ' min')




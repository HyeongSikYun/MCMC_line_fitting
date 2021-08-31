import numpy as np
from MCMC_modules import *
from line_fit_modules import *

log_N_guess = 20
T_ex_guess = 10.0
width_guess = 1.0
Molecule_file = 'e-ch3oh.dat'
T_rms = 0.02
jump_scale = 0.03
n_jump = 30000
freq_range = 1.0e-2 ## in GHz; 1MHz

print 'Read the data'
freq = []; spec = []
with open('test_spectrum.txt','r') as ascii_data:
  for line in ascii_data:
    line = line.strip()
    data = line.split()
    freq = np.append(freq,np.float(data[0]))
    spec = np.append(spec,np.float(data[1]))

lev_data, tr_data = read_mol_file(Molecule_file)

sel_idx = np.zeros(np.shape(freq)[0], dtype=bool)
for i in range(tr_data['N_tr']):
  tr_freq = tr_data['freq'][i]
  sel_idx1 = freq > tr_freq-freq_range
  sel_idx2 = freq < tr_freq+freq_range
  sel_idx12 = np.logical_and(sel_idx1,sel_idx2)
  sel_idx = np.logical_or(sel_idx, sel_idx12)

print 'run MCMC fitting'
initial_p = np.array([log_N_guess,T_ex_guess,width_guess])
prior_p = initial_p
jump_count = 0
pri_likehood = log_likelyhood_const_sig(freq,spec,prior_p,T_rms,lev_data,tr_data,sel_idx)
MC_likehood = pri_likehood
MC_log_N=prior_p[0]; MC_T_ex=prior_p[1]; MC_width=prior_p[2]
while jump_count < n_jump:
  pri_likehood = log_likelyhood_const_sig(freq,spec,prior_p,T_rms,lev_data,tr_data,sel_idx)
  trial_dp = rand_norm(0.0, jump_scale,3)
  trial_p = prior_p + trial_dp
  tri_likehood = log_likelyhood_const_sig(freq,spec,trial_p,T_rms,lev_data,tr_data,sel_idx)
  judge = judge_acceptance(10.0**(tri_likehood-pri_likehood))
  if judge == True:
    jump_count = jump_count + 1
    print "{0:4.1f} %".format(np.float(jump_count)/np.float(n_jump)*100)
    prior_p = trial_p
    MC_likehood = np.append(MC_likehood,tri_likehood)
    MC_log_N = np.append(MC_log_N,trial_p[0])
    MC_T_ex = np.append(MC_T_ex,trial_p[1])
    MC_width = np.append(MC_width,trial_p[2])

print 'generate output file'
out_file = open('out_MC_test.txt','wb')
for i in range(np.shape(MC_likehood)[0]):
  string = "{0:}    {1:}    {2:}    {3:}\n".format(MC_likehood[i],MC_log_N[i],MC_T_ex[i],MC_width[i])
  out_file.write(string)
out_file.close()

print 'Finished'
print 'use the other code to check the output file'

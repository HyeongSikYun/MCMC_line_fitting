import numpy as np
import matplotlib.pyplot as plt

def read_mol_file(name):
  h = 6.62607004E-034
  c = 299790000
  k = 1.38064852E-023
  with open(name,'r') as ascii_data:
    count = 0
    for i in range(5):
      line = ascii_data.readline()
    line = ascii_data.readline()
    line = line.strip()
    data = line.split()
    n_levs = int(data[0])
    line = ascii_data.readline()
    lev_idx = []; lev_E_K = []; lev_weight = []; lev_quantum = []
    for i in range(n_levs):
      line = ascii_data.readline()
      line = line.strip()
      data = line.split()
      lev_idx = np.append(lev_idx,np.int(data[0]))
      lev_E_K = np.append(lev_E_K,np.float(data[1])*100.0*h*c/k)
      lev_weight = np.append(lev_weight,np.float(data[2]))
      lev_quantum = np.append(lev_quantum,data[3])
    print str(n_levs)+' levels are readed'
    lev_data = {}
    lev_data['N_lev'] = n_levs
    lev_data['idx'] = lev_idx
    lev_data['E_K'] = lev_E_K
    lev_data['weight'] = lev_weight
    lev_data['quantum_number'] = lev_quantum
    line = ascii_data.readline()
    line = ascii_data.readline()
    line = line.strip()
    data = line.split()
    n_trans = int(data[0])
    line = ascii_data.readline()
    tr_idx = []; tr_up = []; tr_low = []; tr_Ein_A = []; tr_freq = []; tr_E_up = []
    for i in range(n_trans):
      line = ascii_data.readline()
      line = line.strip()
      data = line.split()
      tr_idx = np.append(tr_idx,int(data[0]))
      tr_up = np.append(tr_up,int(data[1]))
      tr_low = np.append(tr_low,int(data[2]))
      tr_Ein_A = np.append(tr_Ein_A,np.float(data[3]))
      tr_freq = np.append(tr_freq,np.float(data[4]))
      tr_E_up = np.append(tr_E_up,np.float(data[5]))
    print str(n_trans)+' transitions are readed'
    tr_data = {}
    tr_data['N_tr'] = n_trans
    tr_data['idx'] = tr_idx
    tr_data['lev_up'] = tr_up
    tr_data['lev_low'] = tr_low
    tr_data['Ein_A'] = tr_Ein_A
    tr_data['freq'] = tr_freq
    tr_data['E_up'] = tr_E_up
  return lev_data, tr_data

def Boltzmann_distribution(Eu,El,gu,gl,T):
  k = 1.38064852E-023
  if Eu <= El:
    print 'Please input the required parageters as follows'
    print 'E_up, E_low, g_up, g_low, T_ex'
    aaa=aaa
  else:
    Eu = Eu*k; El = El*k
    exponent = (El-Eu)/(k*T)
    ratio = (np.float(gu)/np.float(gl)) * np.exp(exponent)
  return ratio

def Thermalized_level_population(Lev_data,T):
  Lev_population = np.zeros(Lev_data['N_lev'])
  Lev_population[0] = 1
  E_l = Lev_data['E_K'][0]
  g_l = Lev_data['weight'][0]
  for i in range(Lev_data['N_lev']-1):
    E_u = Lev_data['E_K'][i+1]
    g_u = Lev_data['weight'][i+1]
    ratio = Boltzmann_distribution(E_u,E_l,g_u,g_l,T)
    Lev_population[i+1] = ratio
  return Lev_population

def parti_func_thermalize(Lev_data,T):
  population = Thermalized_level_population(Lev_data,T)
  Total = np.sum(population)
  return population[0]/Total

def Gaussian(x,cen,sig,height):
  exponent = 0-((x-cen)**2.0/(2*sig**2.0))
  y = height*np.exp(exponent)
  return y

def derive_alpha(T_ex):
  return 0-(np.log10(2.71828182846)/T_ex)

def derive_beta(lev_data,T_ex, N_mol):
  Q = parti_func_thermalize(lev_data,T_ex)
  return np.log10(N_mol/Q)

def make_spectrum(freq, lev_data, tr_data, p):
  h = 6.62607004E-034
  c = 299790000
  k = 1.38064852E-023
  freq = freq*1.0e+9
  N_mol = 10.0**p[0]
  T_ex = p[1]
  width = p[2]
  alpha = derive_alpha(T_ex)
  beta = derive_beta(lev_data,T_ex,N_mol)
### select the transitions within the frequency range ###
  min_freq = np.min(freq); max_freq = np.max(freq)
  line_freqs = tr_data['freq']*1.0e+9
  sel_idx1 = line_freqs > min_freq
  sel_idx2 = line_freqs < max_freq
  sel_idx = np.logical_and(sel_idx1,sel_idx2)
  n_trns = np.sum(sel_idx)
#  print 'among the transitions in molecular data,'
#  print str(n_trns)+' transitions are selected.'
  tr_E_up = tr_data['E_up'][sel_idx]
  tr_lev_up = tr_data['lev_up'][sel_idx]
  tr_Ein_A = tr_data['Ein_A'][sel_idx]
  tr_freq = tr_data['freq'][sel_idx]

### derive the intensity from the equation of rotational diagram ###
  inte_space = np.zeros(np.shape(freq)[0])
  x_vals = tr_E_up
  y_vals = (x_vals * alpha) + beta
  N_mols = np.zeros(n_trns)

  for i in range(n_trns):
    tr_idx = tr_lev_up[i]
    g_idx = lev_data['idx'] == tr_idx
    g_up = lev_data['weight'][g_idx]
    N_mols[i] = (10.0**y_vals[i])*g_up
    
  ints = N_mols * tr_Ein_A * tr_freq * 1.0e+9 * h
#  Gauss_cent = []; Gauss_sig = []; Gauss_hight=[]
  for i in range(n_trns):
    sig = (width*1.0e+3*tr_freq[i]*1.0e+9)/(c*2.3548200450309493)
    cen = tr_freq[i]*1.0e+9
    height = ints[i]/(sig*np.sqrt(2*np.pi))
#    Gauss_cent = np.append(Gauss_cent,cen)
#    Gauss_sig = np.append(Gauss_sig,sig)
#    Gauss_height = np.append(Gauss_height,height)
    y_line = Gaussian(freq,cen,sig,height)
    inte_space = inte_space + y_line

### convert from the intensity to the temperature ###

  temp_space = (inte_space * c**2.0)/(2*k*freq**2.0)
  return temp_space

def line_fit_chisq(x,y,N,T,W,yerr,lev_data,tr_data,sel_idx):
  p = np.array([np.log10(N),T,W])
  fit_y = make_spectrum(x, lev_data, tr_data, p)
  residual = y[sel_idx]-fit_y[sel_idx]
  return np.sum(residual**2.0)/(yerr**2.0)

def log_likelyhood_const_sig(x,y,p,yerr,lev_data,tr_data,sel_idx):
  N = 10.0**p[0]
  T = p[1]
  W = p[2]
#  n_elem = np.sum(sel_idx)
#  c = 0.5*n_elem*(0-np.log10(2*np.pi)-np.log10(yerr**2.0))
  chisq_val = line_fit_chisq(x,y,N,T,W,yerr,lev_data,tr_data,sel_idx)
  return 0-(0.5*chisq_val)

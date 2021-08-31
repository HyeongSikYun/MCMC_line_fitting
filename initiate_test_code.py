import numpy as np
import matplotlib.pyplot as plt
from MCMC_modules import *
from line_fit_modules import *
width = 5.0 ## km/s
T_ex = 30 ## in K
N_mol = 1.0e+17 ## in per beam (?)
Molecule_file = 'e-ch3oh.dat' # The name of molecule datat file
h = 6.62607004E-034
c = 299790000
k = 1.38064852E-023
T_rms = 0.03

p = np.array([np.log10(N_mol),T_ex,width])
freq_space = np.linspace(100.5,101.5,2000)
lev_data, tr_data = read_mol_file(Molecule_file)

### Derive the linear equation for the rotational diagram ###
### The equation would be log y = beta + alpha * E_up_K   ###

alpha = derive_alpha(T_ex)
beta = derive_beta(lev_data,T_ex, N_mol)

### derive the intensity from the equation of rotational diagram ###

temp_space = make_spectrum(freq_space, lev_data, tr_data, p)

### add artificial noise signal to the spectrum ###
noise_signals = rand_norm(0.0, T_rms,np.shape(freq_space)[0])
spectrum = temp_space + noise_signals
#plt.plot(freq_space, spectrum)
#plt.show(); aaa=aaa

### save the data ###
with open('test_spectrum.txt','w') as out_file:
  for i in range(np.shape(freq_space)[0]-1):
    out_file.write('{0:11.5e}    {1:11.5e}'.format(freq_space[i+1],spectrum[i+1])+"\n")

out_file.close()
print 'Finished'

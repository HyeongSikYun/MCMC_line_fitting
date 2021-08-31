import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from MCMC_modules import *

out_file = 'out_MC_test.txt'
Burn_in_steps = 970
show_step_range = [0,Burn_in_steps*2]

likehood=[]; log_N=[]; T_ex=[] ; width = []
with open('./'+out_file,'r') as ascii_data:
  for line in ascii_data:
    line = line.strip()
    data = line.split()
    likehood = np.append(likehood,np.float(data[0]))
    log_N = np.append(log_N,np.float(data[1]))
    T_ex = np.append(T_ex,np.float(data[2]))
    width = np.append(width,np.float(data[3]))

step_idx = np.arange(np.shape(likehood)[0])+1
idx1 = step_idx > show_step_range[0]
idx2 = step_idx < show_step_range[1]
idx = np.logical_and(idx1,idx2)
burn_idx1 = step_idx < Burn_in_steps
burn_idx = np.logical_and(idx,burn_idx1)
post_idx = np.logical_and(idx,~burn_idx1)

fig = plt.figure(101, figsize=[6.694,6.694])
layout = gridspec.GridSpec(2,2,left=0.11,right=0.98,bottom=0.1,top=0.96,
     wspace=0.25)

first_fig = fig.add_subplot(layout[0,0])
plt.plot(step_idx[burn_idx],likehood[burn_idx],'-',color='red')
plt.plot(step_idx[post_idx],likehood[post_idx],'-',color='blue')
first_fig.set_xlabel('Steps')
first_fig.set_ylabel('log(Likelyhood)')
first_fig.set_xlim(show_step_range)

second_fig = fig.add_subplot(layout[0,1])
plt.plot(log_N[0],T_ex[0],'x ',color='black',markersize=7)
plt.plot(log_N[burn_idx],T_ex[burn_idx],'-',color='red')
plt.plot(log_N[post_idx],T_ex[post_idx],'-',color='blue')
second_fig.set_xlabel('$Log$($N_\mathrm{mol}$)')
second_fig.set_ylabel('T$_\mathrm{ex}$ (K)')

third_fig = fig.add_subplot(layout[1,0])
plt.plot(T_ex[0],width[0],'x ',color='black',markersize=7)
plt.plot(T_ex[burn_idx],width[burn_idx],'-',color='red')
plt.plot(T_ex[post_idx],width[post_idx],'-',color='blue')
third_fig.set_xlabel('T$_\mathrm{ex}$ (K)')
third_fig.set_ylabel('Line width (km s$^{-1}$)')

fourth_fig = fig.add_subplot(layout[1,1])
plt.plot(log_N[0],width[0],'x ',color='black',markersize=7)
plt.plot(log_N[burn_idx],width[burn_idx],'-',color='red')
plt.plot(log_N[post_idx],width[post_idx],'-',color='blue')
fourth_fig.set_xlabel('$Log$($N_\mathrm{mol}$)')
fourth_fig.set_ylabel('Line width (km s$^{-1}$)')

fig2 = plt.figure(102, figsize=[6.694,6.694])
layout2 = gridspec.GridSpec(3,3,left=0.11,right=0.98,bottom=0.1,top=0.9,wspace=0.0,hspace=0.0)

triangle_idx = step_idx > Burn_in_steps
post_log_N = log_N[triangle_idx]; post_T_ex = T_ex[triangle_idx]; post_width = width[triangle_idx]

### generate contour maps
log_N_bins, T_ex_bins, first_im_mat = make_triangle_image(post_log_N,post_T_ex)
log_N_vals = 0.5*(log_N_bins[1:]+log_N_bins[0:-1])
T_ex_vals = 0.5*(T_ex_bins[1:]+T_ex_bins[0:-1])
log_N_mat, T_ex_mat = np.meshgrid(log_N_vals,T_ex_vals)
first_fig = fig2.add_subplot(layout2[1,0])
plt.contour(log_N_mat,T_ex_mat,first_im_mat)
first_fig.xaxis.set_visible(False)

log_N_bins, width_bins, second_im_mat = make_triangle_image(post_log_N,post_width)
log_N_vals = 0.5*(log_N_bins[1:]+log_N_bins[0:-1])
width_vals = 0.5*(width_bins[1:]+width_bins[0:-1])
log_N_mat, width_mat = np.meshgrid(log_N_vals,width_vals)
second_fig = fig2.add_subplot(layout2[2,0])
plt.contour(log_N_mat, width_mat, second_im_mat)

T_ex_bins, width_bins, third_im_mat = make_triangle_image(post_T_ex,post_width)
T_ex_vals = 0.5*(T_ex_bins[1:]+T_ex_bins[0:-1])
width_vals = 0.5*(width_bins[1:]+width_bins[0:-1])
T_ex_mat, width_mat = np.meshgrid(T_ex_vals, width_vals)
third_fig = fig2.add_subplot(layout2[2,1])
plt.contour(T_ex_mat, width_mat, third_im_mat)
third_fig.yaxis.set_visible(False)

###generate histograms
fifth_fig = fig2.add_subplot(layout2[0,0])
plt.hist(post_log_N,log_N_bins,histtype='bar')
mean_val = np.mean(post_log_N)
plt.axvline(x=mean_val,linestyle='--',color='black')
mean_txt = "{0:4.1f}".format(mean_val)
plt.text(0.95,0.95,mean_txt,transform=fifth_fig.transAxes,ha='right',va='top')
fifth_fig.set_xlim([np.min(log_N_vals),np.max(log_N_vals)])
fifth_fig.set_title('$Log$($N_\mathrm{mol}$)')
fifth_fig.xaxis.set_visible(False)
fifth_fig.yaxis.set_visible(False)

sixth_fig = fig2.add_subplot(layout2[1,1])
plt.hist(post_T_ex,T_ex_bins,histtype='bar')
mean_val = np.mean(post_T_ex)
plt.axvline(x=mean_val,linestyle='--',color='black')
mean_txt = "{0:4.1f}".format(mean_val)
plt.text(0.95,0.95,mean_txt,transform=sixth_fig.transAxes,ha='right',va='top')
sixth_fig.set_xlim([np.min(T_ex_vals),np.max(T_ex_vals)])
sixth_fig.set_title('T$_\mathrm{ex}$ (K)')
sixth_fig.xaxis.set_visible(False)
sixth_fig.yaxis.set_visible(False)

seventh_fig = fig2.add_subplot(layout2[2,2])
plt.hist(post_width,width_bins,histtype='bar')
mean_val = np.mean(post_width)
plt.axvline(x=mean_val,linestyle='--',color='black')
mean_txt = "{0:4.1f}".format(mean_val)
plt.text(0.95,0.95,mean_txt,transform=seventh_fig.transAxes,ha='right',va='top')
seventh_fig.set_xlim([np.min(width_vals),np.max(width_vals)])
seventh_fig.set_title('Line width (km s$^{-1}$)')
seventh_fig.yaxis.set_visible(False)

plt.show(); aaa=aaa

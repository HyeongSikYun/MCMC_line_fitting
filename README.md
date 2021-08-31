# MCMC_line_fitting
This is an application of the MCMC linear fitting code to derive the gas column density and excitation temperature from the observed molecular line spectra. 

This code is based on the Markov Chain Monte Carlo (MCMC) class in Sagan Summer Workshop in 2016 presented by David Kipping. The class video is accessible via an youtube channel, 'Sagan Summer Workshop' Please use the following link https://www.youtube.com/watch?v=vTUwEu53uzs

This is my first try to use the MCMC procedure in an astronomical subject.
It contains several python codes as follows: (1) a code to generate artificial line spectra for a given molecule, (2) a code to run the rotational diagram analysis using MCMC procedure, and (3) draw the triangle plot of the MCMC result. you can obtain the molecular data file from Leiden Atomic and Molecular Database (LAMDA; https://home.strw.leidenuniv.nl/~moldata/). 

### How to use
### 1. Open 'initiate_test_code.py' and set the true values.
### 2. run initiate_test_code.py
### 3. run test_MCMC_run.py
### 4. run draw_triangle_plot.py

Followings are the description of the files.
>> MCMC_modules.py: it contains several subroutines to generate random number and to judge the acceptance of a randomized step.
>> line_fit_modules.py: it contains several subroutines to run the MCMC test for the rotational diagram method.
>> initiate_test_code.py: the code to generate the artificial spectrum to test the rotational diagram method with the MCMC procedure.
>> test_spectrum.txt: the output of 'initiate_test_code.py', which contains the artificial spectrum
>> test_MCMC_run.py: the code to run the MCMC procedure to find the column density and excitation temperature of line emitting molecules.
>> out_MC_test.txt: the output file of the MCMC procedure
>> draw_triangle_plot.py: the code to generate a triangle plot to check the result of the MCMC code. You would change the 'show_step_range' and 'Burn_in_steps'.

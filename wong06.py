#!/usr/bin/env python
"""
Simulation for two-populatio reduced mean-field attractor model
(Wong K-F & Wang X-J, J Neurosci (2006))
Author: John D. Murray (john.murray@aya.yale.edu)
Updated: 2013-07-13
"""

import numpy as np
import matplotlib.pyplot as plt

### Stimulus parameters

coh         = 100. # Motion coherence (positive sign is toward RF of pop1) [%]
mu0         = 30.0 # Zero motion coherence stimulus strength
noise_amp   = 0.02 # Noise amplitude into selective populations
N_trials    = 1 # Total number of trials
Tstim       = 1000 # Time of stimulus onset [ms]
Tdur        = 1000 # Duration of stimulus [ms]

### Network parameters

JN11 = 0.2609 # self-excitation
JN22 = 0.2609 # self-excitation
JN12 = -0.0497 # cross-inhibition
JN21 = -0.0497 # cross-inhibition

Ib1 = 0.3255 # background current
Ib2 = 0.3255 # background current

### Synaptic time and other constants

Tnmda       = 100   # NMDA decay time constant
Tampa       = 2      # AMPA decay time constant
gamma       = 0.641 # Parameter that relates presynaptic input firing rate to synaptic gating variable
JAext       = 0.00052 # Synaptic coupling constant to external inputs

### F-I curve parameters
a = 270
b = 108
d = 0.1540

### Time conditions
dt         = 0.5 # Time step (ms)
time_wind  = 5/dt # Temporal window size for averaging
T_total    = int(5000/dt + time_wind) # Total number of steps

slide_wind = 20/dt # Sliding step for window

# Make lists to store firing rate (r) and gating variable (s)
r1_traj = [];  r2_traj = []; r3_traj = []; rI_traj = []
s1_traj = [];  s2_traj = []; s3_traj = []; sI_traj = []

for i in xrange(N_trials): #Loop through trials

    print "trial # ", i

    #Set random seed
    np.random.seed(i)

    #Initialize
    nu1_wind = [] ; nu2_wind = [] ;
    s1_wind  = [] ; s2_wind  = [] ;

    cross_r1 = []; cross_r2 = [];

    s1_in=0.1
    s2_in=0.1
    nu1_in = 2
    nu2_in = 2
    I_eta1_in = noise_amp*np.random.randn(T_total)
    I_eta2_in = noise_amp*np.random.randn(T_total)

    s1 = s1_in*np.ones(T_total)
    s2 = s2_in*np.ones(T_total)
    nu1 = nu1_in*np.ones(T_total)
    nu2 = nu2_in*np.ones(T_total)
    phi1 = nu1_in*np.ones(T_total)
    phi2 = nu2_in*np.ones(T_total)
    I_eta1 = I_eta1_in*np.ones(T_total)
    I_eta2 = I_eta2_in*np.ones(T_total)

    Isyn1 = np.zeros(T_total)
    Isyn2 = np.zeros(T_total)

    for t in xrange(0,T_total-1): #Loop through time for a trial

        #---- Random dot stimulus------------------------------------------------------

        I_stim_1 = (Tstim/dt<t & t<(Tstim+Tdur)/dt)*(JAext*mu0*(1+coh/100)); # To population 1
        I_stim_2 = (Tstim/dt<t & t<(Tstim+Tdur)/dt)*(JAext*mu0*(1-coh/100)); # To population 2

        #---- Resonse function of competiting excitatory population 1 ------

        # Total synaptic input to population 1
        Isyn1[t] = JN11*s1[t] + JN12*s2[t] + I_stim_1 + I_eta1[t];

        # Transfer function of population 1
        phi1[t]  = (a*Isyn1[t]-b)/(1-np.exp(-d*(a*Isyn1[t]-b)));

        #---- Response function of competiting excitatory population 2 -----

        # Total synaptic input to population 2
        Isyn2[t] = JN22*s2[t] + JN21*s1[t] + I_stim_2 + I_eta2[t];

        # Transfer function of population 2
        phi2[t]  = (a*Isyn2[t]-b)/(1-np.exp(-d*(a*Isyn2[t]-b)));

        #---- Dynamical equations -------------------------------------------

        # Mean NMDA-mediated synaptic dynamics updating
        s1[t+1] = s1[t] + dt*(-s1[t]/Tnmda + (1-s1[t])*gamma*nu1[t]/1000);
        s2[t+1] = s2[t] + dt*(-s2[t]/Tnmda + (1-s2[t])*gamma*nu2[t]/1000);

        # Ornstein-Uhlenbeck generation of noise in pop1 and 2
        I_eta1[t+1] = I_eta1[t] + (dt/Tampa)*(Ib1-I_eta1[t]) + np.sqrt(dt/Tampa)*noise_amp*np.random.randn() ;
        I_eta2[t+1] = I_eta2[t] + (dt/Tampa)*(Ib2-I_eta2[t]) + np.sqrt(dt/Tampa)*noise_amp*np.random.randn() ;

        # To ensure firing rates are always positive. Large noise amplitude
        # may result in unwanted negative values
        if phi1[t] < 0:
            nu1[t+1] = 0
            phi1[t] = 0
        else:
            nu1[t+1] = phi1[t]
        if phi2[t] < 0:
            nu2[t+1] = 0
            phi2[t] = 0
        else:
            nu2[t+1] = phi2[t]

    # ---- Calculating the mean rates and gating variables with sliding window -----

    nu1_wind.append(np.mean(nu1[1:time_wind]))
    nu2_wind.append(np.mean(nu2[1:time_wind]))

    s1_wind.append(np.mean(s1[1:time_wind]))
    s2_wind.append(np.mean(s2[1:time_wind]))

    for j in xrange(int((T_total-time_wind)/slide_wind)):
        nu1_wind.append(np.mean(nu1[j*slide_wind:j*slide_wind+time_wind]))
        nu2_wind.append(np.mean(nu2[j*slide_wind:j*slide_wind+time_wind]))
        s1_wind.append(np.mean(s1[j*slide_wind:j*slide_wind+time_wind]))
        s2_wind.append(np.mean(s2[j*slide_wind:j*slide_wind+time_wind]))

    r1_traj.append(nu1_wind)
    r2_traj.append(nu2_wind)
    s1_traj.append(s1_wind)
    s2_traj.append(s2_wind)

### Vectorizing variables for evaluations at end of (block) loop

#Plot
n_plot = 10
plt.figure()
for w in xrange(min(n_plot,N_trials)):
    print 'w', w
    print np.shape(r1_traj)
    print np.shape(np.arange(dt*time_wind-dt*slide_wind,dt*T_total-dt*time_wind,dt*slide_wind))
    plt.plot(np.arange(dt*time_wind-dt*slide_wind,dt*T_total,dt*slide_wind),r1_traj[w],c='b')
    plt.plot(np.arange(dt*time_wind-dt*slide_wind,dt*T_total,dt*slide_wind),r2_traj[w],c='r')
plt.xlim(left=0)
plt.xlabel('Time (ms)')
plt.ylabel('Firing rate (Hz)')
plt.show()
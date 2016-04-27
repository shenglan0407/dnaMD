#! /usr/bin/env python

##############################################################################
# Copyright 2015 Stanford University and the Author
#
# Author: Shenglan Qiao
#
#
#
#############################################################################


##############################################################################
# Imports
##############################################################################

import numpy as np
import time

import h5py
import os

import matplotlib.pyplot as plt
##############################################################################
# Code
##############################################################################

dna_model = "dna17bp_sol_neu-MD_3_sliced"

f = h5py.File("simulated_data/dna17bp_sol_neu-MD_3_sliced_1_10000.hdf5",'a')
shots = set(f.keys())
shots.remove('phi_values')
shots.remove('q_values')

n_shots=len(shots)
print "there are %d exposure in this file"%n_shots
# load PI
PI=[]
for ii in shots:
    PI.append(f[ii][0][0])
phi=f['phi_values'][:]
q_values=f['q_values'][:]

PI = np.array(PI)
ave_PI = np.mean(PI,axis=0)
for ii in range(n_shots):
    PI[ii]=PI[ii]-ave_PI
    
stride=1


def calc_corr(shot, delta):
    
    return np.mean([shot[ii]*shot[ii+delta] for ii in range(len(phi)-delta)])
    
        

deltas = np.arange(0,len(phi),stride)
corr = np.zeros(len(deltas))

tic = time.clock()
for this_shot in PI:
    corr =corr+[calc_corr(this_shot,delta) for delta in deltas]
toc = time.clock()
print "Time for computing corr: %.2f"%(toc-tic)

# save the correlation data in the same hd5f file
f.create_dataset('intercorr',data=corr)
f.close()

#show a plot
# plt.plot(deltas,corr)
# plt.xlim(2,360)
# plt.xlabel('phi(deg)')
# plt.ylabel('Corr')
# plt.show()


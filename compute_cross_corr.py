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

import h5py
import os

import sys
import getopt


import matplotlib.pyplot as plt
##############################################################################
# Code
##############################################################################

def calc_corr(shot1,shot2, delta,n_phi):
    
    return np.mean([shot1[ii]*shot2[ii+delta] for ii in range(n_phi-delta)])

def usage():
    print './compute_cross_corr.py -i <PI> -n <n_files> -o <outputfile> -s <q1> -e <q2>'

def main(argv):
##############################################################################
# default settings
    n_files = 50
    input_series = 'thor_run13_PI'
    output_file = 'simulated_data/50.0_m1_crosscorr.hdf5'
    q1id=51
    q2id=71
##############################################################################

    try:
        opts, args = getopt.getopt(argv,"hi:n:o:s:e:",["PI=","n_files=","ofile=","q1=","q2="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--PI"):
            input_series = arg
        elif opt in ("-n", "--n_files"):
            n_files = int(arg)
        elif opt in ("-o", "--ofile"):
            output_file = "simulated_data/%s.hdf5"%arg
        elif opt in ("-s", "--q1"):
            q1id=int(arg)
        elif opt in ("-e", "--q2"):
            q2id=int(arg)

    



    files = ['simulated_data/%s_%d.hdf5'%(input_series,num) for num in range(n_files)]
    test_data=h5py.File(files[0],'r')
    stride=2
    phi=test_data['phi_values'][:]
    n_phi=len(phi)
    q_values=test_data['q_values'][:]
    deltas = np.arange(0,len(phi),stride)
    corr = np.zeros(len(deltas))

    test_data.close()

##############################################################################
# print settings  
    print "Computing cross-corr between q1 = %.3f and q2 = %.3f"%(q_values[q1id],q_values[q2id])
    print "Output saved in %s\n" %output_file

##############################################################################
    output_data=h5py.File(output_file,'a')
    output_label='%.3f,%.3f'%(q_values[q1id],q_values[q2id])
    if output_label in output_data.keys():
        print "ERROR: cross_corr for this pairs of q values already exist. Please double-check what you are doing!"
        sys.exit()

    n_shots=0
    for this_file in files:
        data=h5py.File(this_file,'r')
    
        PI1=data['polar_intensities'][:,q1id,:]
        PI2=data['polar_intensities'][:,q2id,:]
    
        n_shots=n_shots+len(PI1)
    #     print n_shots
    #     print len(PI1[0])
    
        for shots in zip(PI1,PI2):
            corr =corr+[calc_corr(shots[0],shots[1],delta,n_phi) for delta in deltas]
    
    
        data.close()

#     print n_shots
    corr=np.array(corr)/n_shots 
    
    output_data.create_dataset(output_label,data=corr)
    output_data[output_label].attrs.create('n_shots',n_shots)
    output_data.close()
    
#     plt.plot(deltas,corr)
#     plt.show()



if __name__ == "__main__":
   main(sys.argv[1:])

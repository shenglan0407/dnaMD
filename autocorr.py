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
import time

##############################################################################
# Code
##############################################################################
def usage():
    print './autocorr.py -i <PI_input> -o <outputfile> -a <avePI>'

def calc_corr(shot, delta):
    return np.mean([shot[ii]*shot[ii+delta] for ii in range(len(shot)-delta)])
    
    
def main(argv):
##############################################################################
# default settings
    input_file = None
    output_file = None
##############################################################################


    try:
        opts, args = getopt.getopt(argv,"hi:o:a:",["PI_input=","ofile=","avePI="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--PI_input"):
            input_file = "simulated_data/PI/%s.hdf5"%arg
        elif opt in ("-o", "--ofile"):
            output_file = "simulated_data/autocorr/%s.hdf5"%arg
        elif opt in ("-a", "--avePI"):
            avePI_file = "simulated_data/PI/%s.hdf5"%arg
    
    if input_file == None or avePI_file == None:
        print "Must provide input file and average PI file"
        sys.exit()
        
    if output_file is None:
        output_file ="simulated_data/autocorr/%s"%(input_file.split('PI_')[0].split('PI/')[-1]+'autocorr_'+input_file.split('PI_')[1])
    while os.path.isfile(output_file):
        print "Will not overwrite old file. Please enter new name:"
        output_file = "simulated_data/autocorr/%s.hdf5"%raw_input()
##############################################################################
# print settings  
    print "Computing autocorr for simulated PI %s"%input_file
    print "Output saved in %s\n" %output_file
    

##############################################################################
    
    try:
        f = h5py.File(avePI_file,'r')
        
    except IOError:
        print "average PI file does not exist"
        sys.exit()
    avePI = f['polar_intensities'][:]
    phi=f['phi_values'][:]
    q_values=f['q_values'][:]
    n_q = len(q_values)
    
    f.close()
    
    try:
        f = h5py.File(input_file,'r')
        PI=f['polar_intensities'][:]
        n_shots = PI.shape[0]
#         n_shots = 1

        f.close()
    except IOError:
        print "input file %s does not exist. quit!"%input_file
        sys.exit()
    
    stride = 1
    deltas = np.arange(0,len(phi),stride)
    autocorr = np.zeros((PI.shape[0],n_q,len(deltas)))
    print "compute autocorrelator for %d shots"%n_shots
    print "and %d q values"%n_q
    
    tic = time.clock()
    for idx in range(n_shots):
        PI[idx] = PI[idx]-avePI
        for q_idx in range(n_q):
            autocorr[idx,q_idx,:] = np.array([calc_corr(PI[idx,q_idx,:],delta) for delta in deltas])
    toc = time.clock()
    print "Time for computing autocorr: %.2f"%(toc-tic) 
    
    f = h5py.File(output_file,'w')
    f.create_dataset('autocorr',data=autocorr)
    f.create_dataset('q_values', data=q_values)
    f.create_dataset('phi_values',data=phi[deltas])

    f.close()



if __name__ == "__main__":
   main(sys.argv[1:])
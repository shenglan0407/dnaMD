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
import mdtraj
import h5py
import os
import sys
import getopt
import time

from thor import xray
import matplotlib.pyplot as plt
##############################################################################
# Code
##############################################################################
def usage():
    print './polar_intensities.py -i <mol_model> -o <outputfile> -q <qvalues> -p <nphi> -s <frames>'

def main(argv):
##############################################################################
# default settings
    dna_model = "dna17bp_sol_neu-MD_1_sliced"

    # default simulation parameters
    frames = [0,10]
    n_shots = 1                      # total number of shots to do
    n_molecules = 1                     # the number of molecules to include per shot
    q_values = [1.71]  # the |q| values of the rings to sim
    n_phi = 360                         # number of pts around the rings
    
    output_file = None
##############################################################################


    try:
        opts, args = getopt.getopt(argv,"hi:o:q:p:s:",["mol_model=","ofile=","qvalues=","nphi=","frames=","nmol="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--mol_model"):
            dna_model = arg
        elif opt in ("-o", "--ofile"):
            output_file = "simulated_data/%s.hdf5"%arg
        elif opt in ("-q", "--qvalues"):
            qs = arg.split('/')
            
            if len(qs) > 2:
                try:
                    q_values = np.arange(float(qs[0]),float(qs[1]),float(qs[2]))
                except ValueError:
                    print 'Enter lower bound, upper bound, step for q values, separated by /'
                    sys.exit(2)
            elif len(qs) == 2:
                try:
                    q_values = np.array([float(c) for c in qs])
                except ValueError:
                    print 'Enter two q values, separated by /'
                    sys.exit(2)
            else:
                try:
                    q_values = [float(arg)]
                except ValueError:
                    print 'Enter a single value of q of a range of values separated by space'
                    sys.exit(2)
                
        elif opt in ("-p","--nphi"):
            n_phi = int(arg)
        
        elif opt in ("-s","--frames"):
            frames = arg.split('/')
            
            if len(frames)!=2:
                print "Please enter the starting frame and the ending frame numbers separated by /"
                sys.exit(2)
            else:
                frames = [int(c) for c in frames]
                if frames[1]<frames[0]:
                    print "Starting frame must be smaller than the ending one"
                    sys.exit(2)
                    
    
    t=mdtraj.load("dna_models/%s.trr"%dna_model,top="dna_models/%s.pdb"%dna_model)[frames[0]:frames[1]]
    if output_file is None:
        output_file ="simulated_data/%s_%d_%d.hdf5"%(dna_model,frames[0],frames[1])
    while os.path.isfile(output_file):
        print "Will not overwrite old file. Please enter new name:"
        output_file = "simulated_data/%s.hdf5"%raw_input()
##############################################################################
# print settings  
    print "Computing polar intensities for model %s"%dna_model
    print "Here is some info about the trajactories"
    print t
    print "Size of the simulation box: %.3f %.3f %.3f"%(t.unitcell_lengths[0][0],t.unitcell_lengths[0][1],t.unitcell_lengths[0][2])
    print "--> Frames to compute ranges from %d to %d"%(frames[0],frames[1])
    print "--> q values: %s"%str(q_values)
    print "--> Num of phi between 0 and pi: %d"%n_phi
    print "Output saved in %s\n" %output_file
    

##############################################################################
# do the simulations
    f = h5py.File(output_file,'a')
#     idx = frames[0]
    idx = 0
    tic = time.clock()
    all_PI = np.zeros((t.n_frames,len(q_values),n_phi))
    for this_frame in t:

#         print this_frame
        rings = xray.Rings.simulate(this_frame, n_molecules, q_values, n_phi, n_shots)

        # data to save
        PI=rings.polar_intensities
        
        # if it is the first frame, save phi, q data
        if idx == 0:
            phis=rings.phi_values
            qs=rings.q_values
            
            f.create_dataset('phi_values',data=phis)
            f.create_dataset('q_values',data=qs)

        # save the data
        all_PI[idx]=PI
        idx+=1
    # save the data
    f.create_dataset('polar_intensities',data=all_PI)
    f.close()
    
    toc=time.clock()
    print "Total computation time is %.3g"%(toc-tic)


if __name__ == "__main__":
   main(sys.argv[1:])

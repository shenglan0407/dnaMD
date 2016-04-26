import numpy as np
import os

ice_model = "50.0_m1"
n_runs = 10 # runs per batch
n_batches = 5
n_shots = 200
q_low = 1.35
q_up = 2.15
step = 0.005


for this_batch in range(n_batches):
    with open(os.getcwd()+'/barley_scripts/submitBatch%d.sh'%this_batch,'w') as submit:
        for this_run in range(n_runs*n_batches)[this_batch*n_runs:(this_batch+1)*n_runs]:
            sub_file = 'job_%d.sh'%this_run
            output = 'thor_run13_PI_%d'%this_run
            with open(os.getcwd()+'/barley_scripts/'+sub_file,'w') as this_file:
                this_file.write('#!/bin/bash\n')
                this_file.write('#$ -N PI_%d\n'%this_run)
                this_file.write('#$ -wd /srv/zfs01/user_data/shenglan/ice_nanocrystals/\n')
                this_file.write('#$ -M shenglan@stanford.edu\n')
                this_file.write('#$ -m besa\n')
                this_file.write('#$ -j y\n')
                this_file.write('#$ -o /srv/zfs01/user_data/shenglan/ice_nanocrystals/barley_scripts/out_'+str(this_run)+'.log\n')
                this_file.write('python polar_intensities.py -i %s -o %s -q %.3f/%.3f/%.3f -p 360 -s %d -n 1\n'\
                %(ice_model,output,q_low,q_up,step,n_shots))
                this_file.write('hostname\n')
            submit.write('qsub '+sub_file+'\n')
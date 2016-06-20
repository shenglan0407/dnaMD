import numpy as np
import os

dna_model = "dna17bp_sol_neu-MD_2_sliced"
n_subs = 50 # num of jobs to run
n_shots = 200
q_low = 0.10
q_up = 0.70
step = 0.01



with open(os.getcwd()+'/barley_scripts/submitPIBatch.sh','w') as submit:
    for this_run in range(n_subs):
        sub_file = 'job_2_%d.sh'%this_run
        output = '%s_PI_%d'%(dna_model,this_run)
        frame_start = this_run*n_shots
        with open(os.getcwd()+'/barley_scripts/'+sub_file,'w') as this_file:
            this_file.write('#!/bin/bash\n')
            this_file.write('#$ -N PI_%d\n'%this_run)
            this_file.write('#$ -wd /srv/zfs01/user_data/shenglan/dnaMD/\n')
            this_file.write('#$ -M shenglan@stanford.edu\n')
            this_file.write('#$ -m besan\n')
            this_file.write('#$ -j y\n')
            this_file.write('#$ -o /srv/zfs01/user_data/shenglan/dnaMD/barley_scripts/out_2_'+str(this_run)+'.log\n')
            this_file.write('python polar_intensities.py -i %s -o %s -q %.3f/%.3f/%.3f -p 360 -s %d/%d \n'\
            %(dna_model,output,q_low,q_up,step,frame_start,frame_start+n_shots))
            this_file.write('hostname\n')
        submit.write('qsub '+sub_file+'\n')

import numpy as np
import os
# 
# dna_model = "dna17bp_sol_neu-MD_3_sliced"
# n_subs = 50 # num of jobs to run
# n_shots = 200
# q_low = 0.70
# q_up = 2.15
# step = 0.01

file_list = ['dna17bp_sol_neu-MD_2_sliced_PI_%d'%i for i in range(50)]
file_list.extend(['dna17bp_sol_neu-MD_3_sliced_PI_%d'%i for i in range(50)])
run_num =0


with open(os.getcwd()+'/barley_scripts/submitPIBatch.sh','w') as submit:
    for this_PI in file_list:
        sub_file = 'job_autocorr_%d.sh'%run_num
        
        with open(os.getcwd()+'/barley_scripts/'+sub_file,'w') as this_file:
            this_file.write('#!/bin/bash\n')
            this_file.write('#$ -N autocorr_%d\n'%run_num)
            this_file.write('#$ -wd /srv/zfs01/user_data/shenglan/dnaMD/\n')
            this_file.write('#$ -M shenglan@stanford.edu\n')
            this_file.write('#$ -m besan\n')
            this_file.write('#$ -j y\n')
            this_file.write('#$ -o /srv/zfs01/user_data/shenglan/dnaMD/barley_scripts/out_autocorr_%d.log \n'%run_num)
            this_file.write('./autocorr.py -i %s -a dna17bp_sol_neu-MD_23_sliced_avePI \n'\
            %this_PI)
            this_file.write('hostname\n')
        submit.write('qsub '+sub_file+'\n')
        run_num += 1
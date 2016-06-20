import h5py
import numpy as np
import matplotlib.pyplot as plt

#generate file list
dna_model2 = 'dna17bp_sol_neu-MD_2_sliced_PI'
dna_model3 = 'dna17bp_sol_neu-MD_3_sliced_PI'
file_list = ['simulated_data/autocorr/dna17bp_sol_neu-MD_2_sliced_autocorr_%d.hdf5'%i for i in range(50)]
file_list.extend(['simulated_data/autocorr/dna17bp_sol_neu-MD_3_sliced_autocorr_%d.hdf5'%i for i in range(50)])

average_polar = np.zeros((145,360))
total_shots = 0
counter=0
for file in file_list:
    counter+=1
    try:
        data = h5py.File(file,'r')
    except IOError:
        print "$s does not exist"%file
        continue
    
    average_polar+= np.sum(data['autocorr'],axis=0)
    total_shots+=data['autocorr'].shape[0]
    if counter == len(file_list):
        output = h5py.File('simulated_data/autocorr/dna17bp_sol_neu-MD_23_sliced_aveAutocorr.hdf5','w')
        output.create_dataset('autocorr',data=average_polar/total_shots)
        output.create_dataset('q_values', data=data['q_values'])
        output.create_dataset('phi_values',data=data['phi_values'])
        output.close()
    data.close()

# plt.imshow(average_polar)
# plt.savfig('simulated_data/average_polar_intensity.png')

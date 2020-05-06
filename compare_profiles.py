# Comparing profiles in a window around features (TSS, TES, etc.):

import pickle
import numpy as np
import HTSeq
from matplotlib import pyplot as plt
import os

# compare profiles
halfwinwidth = 3000
# data_folder = "/home/eranko/Runs/CMV/CMV2_HA-IE1_output/BAM_files_Bowtie2" # Folder containing the profile pickles
data_folder = 'c:/Users/Eran/Desktop/CMV Temp stuff'
files = ['TSS_positions_3000_profile.p', 'TES_positions_3000_profile.p', 'start_codons_3000_profile.p', 'stop_codons_3000_profile.p']


# fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)



for i, f in enumerate(files):
    datafile = os.path.join(data_folder,f)
    profile = pickle.load(open(datafile, 'rb'))
    coverage_sum = profile.mean()
    norm_profile = profile*1.0/coverage_sum   
    ax = plt.subplot(2,2,i)
    ax.plot(np.arange( -halfwinwidth, halfwinwidth ), norm_profile, 'b')
    plt.ylim([0,2])
    ax.set_title(f[:-2])
    
plt.show()


#####################################################
# from Utils import Write, Load
# Write(filename, data)
# Load(filename)


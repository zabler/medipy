'''SCRIPT ECOG1 FREQUENCY'''

'''SYSTEM AND IO'''
import argparse
import os

'''SCIPY AND ANALYSIS '''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''DATA LOCATION'''

# # COMMAND LINE AS GUI
# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--input', dest='input', required=True, type=str, help='Path to input folder')
# args = parser.parse_args() #args ist Objekt mit den gegebenen Informationen aus CMD
# root = args.input

# DIRECT PATH FOR DEBUGGING
root ='C:\\Users\\nzabler\\Desktop\\Testdaten\\Software\\Milli\\03_CSV\\2019_11_29\\oneECOGandpower\\2019-11-29-14-13-30_bic_data_oneECOGandpower.csv'

'''READ'''
# Channel Parser
# Root is a bic file
data = np.genfromtxt(root, delimiter=',',usecols=(1))

'''CONVERSION'''
ecog1 = -0.073 * (data - 3000000 / (2.0 * 46.00))

'''CLEAR'''
#FUTURE Do alway odd, zero padded, windowed,... own fft=?

'''ANALYSIS'''
ecog1_inter = pd.Series(ecog1).interpolate(method='nearest')
N = len(ecog1_inter)
fs = 1000 
spektalraum = 2*abs(np.fft.fft(ecog1_inter))/N
ecog_spec_y = spektalraum[:round(N/2)]
ecog_spec_x = np.linspace(0, fs/2, round(N/2), endpoint=True)

'''PLOTS'''
fig = plt.figure()
fig.suptitle('ECoG')

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(ecog1, lw=0.2,color='black')
ax1.set_title('Time Domain')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Channel1 [µV]')

ax2.plot(ecog_spec_x,ecog_spec_y,lw=0.7,color='black')
ax2.set_title('Frequency Domain')
ax2.set_xlabel('Frequency [Hz]')
ax2.set_ylabel('Magnitude [Energie]')

#plt.savefig(debugpath + '_ecog1.png', dpi=600, bbox_inches='tight')
plt.grid(True)
plt.show()


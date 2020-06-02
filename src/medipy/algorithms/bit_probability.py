'''SCRIPT BIT PROBABILITY'''

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
root ='C:\\Users\\nzabler\\Desktop\\Testdaten\\Software\\Milli\\03_CSV\\2019_11_29\\measurementdata\\2019-11-29-14-13-30_bic_data_measurementdata.csv'

'''READ'''

# Channel Parser
# Root is a bic file
data = np.genfromtxt(root,usecols=(6),delimiter=',')

'''CLEAR DATA'''

# Remove NaN
data = data[~np.isnan(data)]

# Numpy float64 to int
data = data.astype(int)

'''FILTER'''
# EVTL HP Filtern damit DC Offset entfernt?

'''ANALYSIS'''
# Calculating Bit Probability

bit = 16
bits = np.zeros(16)
for k in data:
    current = bin(k)[2:].zfill(bit)
    for ind, val in enumerate(current):
        if val == '1':
            bits[ind]+=1
prob_bits = bits/len(data)*100
 

'''PLOTS'''
bit_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
fig1, ax1 = plt.subplots()
ax1.bar(bit_numbers, height=prob_bits, width=0.6, label='High/1')
#ax1.bar(bit_numbers, height=100-prob_bits, bottom=prob_bits,width=0.6, label='Low/0')
ax1.set_title('Bit Probability')
ax1.set_ylabel('Probabilitiy [%]')
ax1.set_xlabel('Bits')
ax1.grid(True)
ax1.set_ylim(0,100)
ax1.set_xticks(np.arange(1,17))
ax1.set_xticklabels(bit_numbers)
ax1.set_xlim(0,17)
ax1.legend(loc='upper right')
fig1.tight_layout()


#plt.savefig(root+'_bit_probability.png', dpi=1300, bbox_inches='tight')
plt.show()



'''SCRIPT ECOG1 PACKETLOSS'''

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
# # Channel Parser
# # Root is a bic file
data = np.genfromtxt(root, delimiter=',',usecols=(1))

'''CONVERSION'''
ecog1 = -0.073 * (data - 3000000 / (2.0 * 46.00))

'''ANALYSIS'''
trueindex = range(0,len(ecog1))
ecog1_inter = pd.Series(ecog1).interpolate(method='nearest')
paketloss = []
for ind_b, val_b in enumerate(ecog1):
    if np.isnan(val_b):
        paketloss.append(ind_b)

'''PLOTS'''
plt.plot(ecog1, markersize=1, lw=0.2, label=f'Test XY',color='black')
plt.plot(paketloss,ecog1_inter.values[paketloss],marker='x',color='red',markersize=0.7,linestyle = 'None',label='Packetloss')
plt.title('ECoG')
plt.ylabel('Channel1 [µV]')
plt.xlabel('Time [sec]')
plt.xlim(0,250000)
xloc,lab = plt.xticks()
plt.xticks(xloc, [0,50,100,150,200,250])
# plt.xlim(0,600000)
# xloc,lab = plt.xticks()
# plt.xticks(xloc, [0,100,200,300,400,500,600])
plt.ylim(-1000,1000)
plt.grid(True)
plt.legend(loc='lower right')
#plt.savefig(debugpath + '_ecog1.png', dpi=600, bbox_inches='tight')
plt.show()


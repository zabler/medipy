'''SCRIPT TRANSIENT RESPONSE'''

'''SYSTEM AND IO'''
import argparse
import os

'''SCIPY AND ANALYSIS '''
import numpy as np
import pandas as pd
import scipy as sc
import matplotlib.pyplot as plt
import datetime as dt

'''DATA LOCATION'''

# # COMMAND LINE AS GUI
# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--input', dest='input', required=True, type=str, help='Path to input folder')
# args = parser.parse_args() #args ist Objekt mit den gegebenen Informationen aus CMD
# root = args.input

# DIRECT PATH FOR DEBUGGING
root ='C:\\Users\\nzabler\\Desktop\\Testdaten\\Software\\Milli\\03_CSV\\2019_11_29\\implantpower\\2019-11-29-14-32-23_bic_data_implantpower.csv'

'''READ'''
# Test Parser
# Root is a bic file
data = np.genfromtxt(root, delimiter=',',usecols=(0,1,2,3,6))
timestamps= data[:,0]
source =  data[:,1] #USE ONLY DATA OF EXTERNAL UNIT =0
coilVoltage = (data[:,2] - 1173)/72.4
targetVoltage = (data[:,3] - 1173)/72.4
duty_cycle = data[:,4]

'''CLEAR DATA'''
remover=[]
for ind_a, val_a in enumerate(source):
    if val_a == 1:
        remover.append(ind_a)

timestamps = np.delete(timestamps,remover)
coilVoltage = np.delete(coilVoltage,remover)
targetVoltage = np.delete(targetVoltage,remover)
duty_cycle = np.delete(duty_cycle,remover)

trueindex = [0]
for b in range(1,len(timestamps)):
    trueindex.append((timestamps[b]-timestamps[b-1])+trueindex[b-1])


'''PLOTS'''
fig, ax1 = plt.subplots()
ax1.plot(trueindex, coilVoltage, 'black', label='Test 2',lw=0.3)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('coilVoltage in V', color='black')
ax1.tick_params('y', colors='black')
#ax1.set_ylim(7.75,8.25)

ax2 = ax1.twinx()
ax2.plot(trueindex, duty_cycle, marker='x',color='red',markersize=1.2,lw=0.6)
ax2.set_ylabel('DutyCycle in %', color='r')
ax2.tick_params('y', colors='r')
ax2.set_ylim(0,100)

plt.title('Transient response Test 2')
#plt.xlim(0,600)
#plt.ylim(37,45)
plt.xlim(0,5000000)
xloc,lab = plt.xticks()
plt.xticks(xloc, [0,1,2,3,4,5])
plt.grid(True)
#plt.legend(loc='lower right')
plt.savefig(root + '_transientresponse.png', dpi=1300, bbox_inches='tight')
fig.tight_layout()
plt.show()


'''SCRIPT COIL VS DUTY'''

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
#root ='C:\\Users\\nzabler\\Desktop\\Testdaten\\Software\\Milli\\03_CSV\\2019_11_29\\oneECOGandpower\\2019-11-29-14-32-23_bic_data_oneECOGandpower.csv'
root = 'C:\\Users\\nzabler\\Desktop\\2020-03-31-15-44-13_bic_data_oneECOGandpower.csv'

'''READ'''

# Channel Parser
# Root is a csv file
data = np.genfromtxt(root, delimiter=',')#,usecols=(5,6))

'''CLEAR DATA'''
duty_cycle = data[:,5]
coilVoltage = (data[:,6] - 1173)/72.4
trueindex = range(0,len(data))

'''PLOTS'''

fig1, ax1 = plt.subplots()
ax1.plot(trueindex, duty_cycle, marker='x',color='red',markersize=0.3,linestyle='None')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('DutyCycle in %', color='r')
ax1.tick_params('y', colors='r')
ax1.set_ylim(40,52)


# ax2 = ax1.twinx()
# ax2.plot(trueindex, coilVoltage, 'black', label='Test 2',lw=0.3)
# ax2.set_ylabel('coilVoltage in V', color='black')
# ax2.tick_params('y', colors='black')
# ax2.set_ylim(7.75,8.25)


#fig1.suptitle('2020-03-31-15-42-20')
ax1.set_title('2020-03-31-15-44-13-First10')
ax1.set_xlim(0,10000)
#ax1.set_xticks(np.arange(0,100,20))
#ax1.set_xticklabels([0,20,40,60,80,100])
ax1.set_xticklabels([0,2,4,6,8,10])

plt.grid(True)
plt.legend(loc='lower right')

fig1.tight_layout()
plt.savefig('2020-03-31-15-44-13-First10.png', dpi=1300, bbox_inches='tight')
plt.show()


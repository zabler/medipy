'''SCRIPT PACKETLOSS'''

'''SYSTEM AND IO'''
import argparse
import os

'''SCIPY AND ANALYSIS '''
import numpy as np
import pandas as pd
import scipy as sc
import matplotlib.pyplot as plt

'''DATA LOCATION'''

# # COMMAND LINE AS GUI
# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--input', dest='input', required=True, type=str, help='Path to input folder')
# args = parser.parse_args() #args ist Objekt mit den gegebenen Informationen aus CMD
# root = args.input

# DIRECT PATH FOR DEBUGGING
root ='C:\\Users\\nzabler\\Desktop\\Testdaten\\Software\\Milli\\03_CSV\\2019_11_29\\oneECOGandpower'

'''READ, ANALYSIS AND PLOTS'''
# Experiment Parser
# Root is a folder with bic files
pathcontent = os.listdir(root) # lists all Tests
for ind_a, val_a in enumerate(pathcontent):
    data = np.genfromtxt(root +'\\'+ pathcontent[ind_a], delimiter=',',usecols=(1))
    packetloss = [0] # Kumulierte Packetloss Liste
    for b in data:
        if np.isnan(b):
            packetloss.append(packetloss[-1]+1)
        else:
            packetloss.append(packetloss[-1])
    plt.plot(packetloss, markersize=1, lw=0.6, label=f'Test {ind_a+1}: Lost {packetloss[-1]} ({round(100*packetloss[-1]/len(packetloss),2)} %)')

'''PLOTS'''
plt.title('Packet Loss')
plt.ylabel('cumulative sum')
plt.xlabel('Time [s]')
plt.xlim(0,250000)
xloc,lab = plt.xticks()
plt.xticks(xloc, [0,100,200,300,400,500,600])
plt.ylim(0,2500)
plt.grid(True)
plt.legend(loc='upper left')
# plt.savefig(root + '_packetloss.png', dpi=1300, bbox_inches='tight')
plt.show()


'''SCRIPT HUMIDITY'''

'''SYSTEM AND IO'''
import argparse
import os

'''SCIPY AND ANALYSIS '''
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc

'''DATA LOCATION'''

# # COMMAND LINE AS GUI
# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--input', dest='input', required=True, type=str, help='Path to input folder')
# args = parser.parse_args() #args ist Objekt mit den gegebenen Informationen aus CMD
# root = args.input

# DIRECT PATH FOR DEBUGGING
root ='C:\\Users\\nzabler\\Desktop\\Testdaten\\Software\\Milli\\03_CSV\\2019_11_29\\implanthealthevent'

'''READ, CONVERSION AND PLOT'''

# Experiment Parser
# Root is a folder with bic files
pathcontent = os.listdir(root) # lists all Tests
for ind,val in enumerate(pathcontent):
    data = np.genfromtxt(root +'\\'+ pathcontent[ind], delimiter=',',usecols=(2))
    data = 100.0 * data / 65535
    data = sc.medfilt(data, kernel_size=21)
    data = data[:-10]
    plt.plot(data, markersize=1, lw=0.6, label=f'Test {ind+1}')

'''PLOTS'''

plt.title('Humidity')
plt.ylabel('Humidity [%]')
plt.xlabel('Time [s]')
plt.xlim(0,250)
plt.ylim(14,19)
plt.grid(True)
plt.legend(loc='lower right')
#plt.savefig(root + '_humidity.png', dpi=1300, bbox_inches='tight')
plt.show()

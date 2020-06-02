'''SCRIPT HUMIDITY OVER YEAR'''

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


'''DATA LOCATION'''

# # COMMAND LINE AS GUI
# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--input', dest='input', required=True, type=str, help='Path to input folder')
# args = parser.parse_args() #args ist Objekt mit den gegebenen Informationen aus CMD
# root = args.input

# DIRECT PATH FOR DEBUGGING
root ='C:\\Users\\nzabler\\Desktop\\Testdaten\\Software\\Milli\\03_CSV\\' 

'''READ'''
# Series of Experiments Parser
# Root is a folder with folders with bic files 
pathcontent = os.listdir(root) # lists all Experiments
criteria = 'implanthealthevent'
data = []
for a in pathcontent: # check every experiment with 
    if criteria in os.listdir(root+a): # if criteria is ok  
        goalfolder = root+a+'\\'+criteria
        goalcontent = os.listdir(goalfolder)
        exp = []
        for ind_b, val_b in enumerate(goalcontent): # Alles CSV dazu einlesen und Daten speichern
            data_in = np.genfromtxt(goalfolder+'\\'+val_b, delimiter=',',usecols=(2))
            exp.append(100.0 * data_in / 65535)
        data.append(np.concatenate(exp))

'''PLOTS'''
#Boxplot
figbp, axbp = plt.subplots()
axbp.set_title(f'Change of humidity in implant over 1 year: +{(np.median(data[-1])-np.median(data[0])).round(2)}%')
axbp.set_ylabel('Humidity [%]')
axbp.set_xlabel('')
axbp.boxplot(data,showfliers=False)
axbp.grid(True)
axbp.set_ylim(15.0,19.0)
axbp.set_xticks(np.arange(len(pathcontent)+1))
xlabels = pathcontent.copy()
xlabels.insert(0,'')
axbp.set_xticklabels(xlabels, rotation=45)
figbp.tight_layout()
#plt.savefig(rootpath + '_humidyear_boxplot.png', dpi=1300, bbox_inches='tight')


#Median Only
figmed, axmed = plt.subplots()
axmed.set_title(f'Change of humidity in implant over 1 year: +{(np.median(data[-1])-np.median(data[0])).round(2)}%')
axmed.set_ylabel('Humidity [%]')
axmed.set_xlabel('')
med = [] 
for ind, val in enumerate(data):
    med.append(np.median(data[ind]))
d = {'Time':pathcontent,'Vals':med}
df = pd.DataFrame(data=d)
df['Time'] = pd.to_datetime(df['Time'],format='%Y_%m_%d')
#pd.to_datetime(df['Time'],format='%Y_%M_%d')
#df = df.set_index(pd.DatetimeIndex(df['Time']))
axmed.plot(df.Time, df.Vals, marker='x',color='red',markersize=4.0,lw=1.0,label='Median')
axmed.grid(True)
axmed.set_ylim(15.0,19.0)
dstart = dt.datetime(2018,12,19)   
dend = dt.datetime(2019,11,29)
axmed.set_xlim(dstart,dend)
axmed.legend(loc='lower right')
figmed.tight_layout()
#plt.savefig(rootpath + '_humidyear_median.png', dpi=1300, bbox_inches='tight')
plt.show()
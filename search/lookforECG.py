from wda.movisens2python import m2pclass
from wda.rdetections import dectclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
lookforECG
'''

# Input: Signalart    
channelname = 'ecg'

# Movisensobjekt einlesen mit geünschter Signalart 
movisensobject = m2pclass.m2pconverter(
    channelname, 'seizures', showtree=True)
# Signalparameter wählen und Signalwerte berechnen
channel = movisensobject.getentry(channelname)
channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
#channel.signal = channel.signal * float(channel.lsbValue)
fs = int(channel.sampleRate)

#Seizures
seizures = movisensobject.getentry('seizures').event
eR= int(movisensobject.getentry('seizures').sampleRate)

#Frequenzverhältnis
ver = fs/eR

# Bereiche um Seizures ausschneiden
size = 0.5*eR

for index, anfall in enumerate(seizures):
    fig = plt.figure(anfall)
    plt.plot(channel.signal[int(anfall-(size*ver)):int(anfall+(2*size*ver))],label=f'Seizure #{index} Index {seizures[index]}')
    plt.axvline(x=size*ver,color='r',linestyle='--')
    plt.plot(size*ver, 0, 'r--', label='Onset')
    plt.legend()
    
plt.show()

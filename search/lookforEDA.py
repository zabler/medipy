from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
lookforECG

(1) Erstellung eines Movisens Objekts mit bestimmtem Signal

(2) Signalwerte mit LsbValue (und Baseline) umrechnen

(3) Bereiche um Seizures berechnen

(4) Signalstücke jeweils einzeln Plotten

'''

# Input: Anzahl der Messungen, Signalart    

channelname = 'eda'

# Movisensobjekt einlesen mit geünschter Signalart 
movisensobject = m2pclass.m2pconverter(
    channelname, 'seizures', showtree=True)

# Signalparameter wählen und Signalwerte berechnen
channel = movisensobject.getentry(channelname)
channel.signal = channel.signal * float(channel.lsbValue)
fs = int(channel.sampleRate)

#Seizures
seizures = movisensobject.getentry('seizures').event
eR= int(movisensobject.getentry('seizures').sampleRate)

#Frequenzverhältnis
ver = fs/eR

# Bereiche um Seizures ausschneiden
size = 5*eR  #+- 5 Sekunden

for index, anfall in enumerate(seizures):
    fig = plt.figure(index)
    plt.plot(channel.signal[int((anfall-size)*ver):int((anfall+3*size)*ver)], label=f'Seizure #{index} Index {seizures[index]}')
    plt.axvline(x=size*ver,color='r',linestyle='--')
    plt.plot(size*ver, 0, 'r--', label='Onset')
    plt.legend()
    
plt.show()

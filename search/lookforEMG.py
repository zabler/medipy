from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
lookforEMG

(1) Erstellung eines Movisens Objekts mit bestimmtem Signal

(2) Signalwerte mit LsbValue (und Baseline) umrechnen

(3) Bereiche um Seizures berechnen

(4) Signalstücke jeweils einzeln Plotten

'''

# Input: Anzahl der Messungen, Signalart    
channelname = 'EMG8'

# Movisensobjekt einlesen mit geünschter Signalart 
movisensobject = m2pclass.m2pconverter(
    channelname, 'bicepsseizures', showtree=True)

# Signalparameter wählen und Signalwerte berechnen
channel = movisensobject.getentry(channelname)
channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
fs = int(channel.sampleRate)

#Seizures
seizures = movisensobject.getentry('bicepsseizures').event
eR= int(movisensobject.getentry('bicepsseizures').sampleRate)

#Frequenzverhältnis
ver = fs/eR

#Size
size = 0.5 *eR

# Bereiche um Seizures ausschneiden
for index, anfall in enumerate(seizures):
    fig = plt.figure(anfall)
    plt.plot(channel.signal[int((anfall-size)*ver):int((anfall+size)*ver)], label=f'Seizure #{index} Index {seizures[index]}')
    plt.axvline(x=size*ver,color='r',linestyle='--')
    plt.plot(size*ver, 0, 'r--', label='Onset')
    plt.legend()
    
plt.show()

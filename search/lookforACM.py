from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

'''
lookforACM

'''

# Input: Signalart    
channelname = 'acc'

# For Schleife über Anzahl der Messungen

# Movisensobjekt einlesen mit gewünschter Signalart 
movisensobject = m2pclass.m2pconverter(
    channelname, 'seizures', showtree=True)

# Signalparameter wählen und Signalwerte berechnen
channel = movisensobject.getentry(channelname)
fs = int(channel.sampleRate)
le = len(channel.signal)

# Subchannels auslesen
accX = channel.signal[0:le:3] * float(channel.lsbValue)
accY = channel.signal[1:le:3] * float(channel.lsbValue)
accZ = channel.signal[2:le:3] * float(channel.lsbValue)

# Marker wählen / Freqverhältnis Samplerate und Eventrate ist bei beiden gleich
seizures = movisensobject.getentry('seizures').event
eR= int(movisensobject.getentry('seizures').sampleRate)

#Frequenzverhältnis
ver = fs/eR

# Bereiche um Seizures ausschneiden
size = 0.25*eR #+- 0.25 Sekunden

for index, anfall in enumerate(seizures):
    fig = plt.figure(anfall)
    plt.plot(0,0,label=f'Seizure #{index} Index {seizures[index]}')
    plt.plot(accX[int(anfall-size*ver):int(anfall+3*size*ver)], label='ACC X', linewidth=0.7, color='black')
    plt.plot(accY[int(anfall-size*ver):int(anfall+3*size*ver)], label='ACC Y', linewidth=0.7, color='blue')
    plt.plot(accZ[int(anfall-size*ver):int(anfall+3*size*ver)], label='ACC Z', linewidth=0.7, color='green')
    plt.axvline(x=size*ver,color='r',linestyle='--')
    plt.plot(size*ver, 0, 'r--', label='Onset')
    plt.legend()
    
plt.show()


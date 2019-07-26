from wda.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt

'''
POSTERPLOT EMG

(1) Erstellung eines Movisens Objekts mit bestimmtem Signal

(2) Signalwerte mit LsbValue und Baseline umrechnen

(3) Bereiche um Seizures berechnen

(4) Signalstücke jeweils einzeln Plotten

'''

# Input: Anzahl der Messungen, Signalart
anzahl_measurements = 1
channelname = 'EMG2'

# For Schleife über Anzahl der Messungen
for k in range(0, anzahl_measurements):
    # Movisensobjekt einlesen mit geünschter Signalart 
    movisensobject = m2pclass.m2pconverter(
        channelname, 'seizures', showtree=True)
    # Signalparameter wählen und Signalwerte berechnen
    channel = movisensobject.getentry(channelname)
    channel.signal = (channel.signal - int(channel.baseline)) * float(channel.lsbValue)
    fs = channel.sampleRate
    seizures = movisensobject.getentry('seizures').event
    # Bereiche um Seizures ausschneiden
    #fig = plt.figure(k)
    for i in seizures:
        fig = plt.figure(i)
        plt.plot(channel.signal[i-2000:i+2000], label=channelname)
        plt.plot(2000, channel.signal[i], 'r+', label='Seizure')
        plt.title(channelname, fontweight="bold")
        plt.xlabel(f'[s] bei Abtastung mit {fs} Hz')
        plt.ylabel(f'Amplitude in {channel.unit}')
        plt.legend()

plt.show()




# Zuerst alle einlesen oder immer Datenbearbeitung vorher?
# Statisches in File oder Dynamisches in Variable Zwischenspeichern der Daten um Seizure
# benötigt erneutes wiedereinlesen, dynamisch benötigt immer den ganzen prozess


# # Plot erstellen
# plt.plot(channel.signal, label='Channel')

# # Marker
# plt.plot(seizures, channel.signal[seizures],
#          'r+', label='Marker')

# # Settings
# plt.title('Title', fontweight="bold")
# plt.xlabel(f'[s] bei Abtastung mit {fs} Hz')
# plt.ylabel(f'Amplitude in {channel.unit}')

# # Bereiche 10-20 Sekunden
# #plt.xlim(20000 / (1 / int(fs)), 30000 / (1 / int(fs)))
# plt.xlim(seizures[2]-1000, seizures[2]+1000)
# plt.legend()
# plt.show()
#und dynamischem Variablenamen#dynamic = {}
    # dynamic[f'movisensobject{k}'] = m2pclass.m2pconverter(
    #     channelname, 'seizures', showtree=True)
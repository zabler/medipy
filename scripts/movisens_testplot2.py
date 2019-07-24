from Module.movisens2python import m2pclass
import numpy as np
from matplotlib import pyplot as plt

'''
EXAMPLE

(1) Erstellung eines Movisens Objekts mit allen Signalen, die verfügbar sind

(2) Signalwerte mit LsbValue und Baseline umrechnen

(2) ECG Signal, Abtastfrequenz, R-Zacken, Marker und Artefakte auslesen

(3) Signal plotten

'''

if __name__ == '__main__':

    # Objekt erstellen mit allen verfügbaren Signalen
    movisensobject = m2pclass.m2pconverter(showtree=True)

    # SignalEntry ECG und ValuesEntry NN_Live auslesen
    channel = movisensobject.getentry('EMG4')
    channel.signal = (channel.signal - int(channel.baseline)) * \
        float(channel.lsbValue)
    fs = channel.sampleRate
    seizures = movisensobject.getentry('seizures').event

    # Plot erstellen
    plt.plot(channel.signal, label='Channel')

    # Marker
    plt.plot(seizures, channel.signal[seizures],
             'b+', label='Marker')

    # Settings
    plt.title('Title', fontweight="bold")
    plt.xlabel(f'[s] bei Abtastung mit {fs} Hz')
    plt.ylabel(f'Amplitude in {channel.unit}')

    # Bereiche 10-20 Sekunden
    #plt.xlim(20000 / (1 / int(fs)), 30000 / (1 / int(fs)))
    plt.xlim(seizures[2]-1000, seizures[2]+1000)
    plt.legend()
    plt.show()

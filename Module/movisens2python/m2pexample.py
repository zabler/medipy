import m2pclass
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

    # Objekt erstellen, mit Signaltyp ECG
    movisensobject = m2pclass.m2pconverter(showtree=True)

    # 'ecg', 'nn_live', 'marker', 'acc', 'artifact', showtree=True)

    # SignalEntry ECG und ValuesEntry NN_Live auslesen
    ecg = movisensobject.getentry('EMG4')
    ecg_signal = (ecg.signal - int(ecg.baseline)) * float(ecg.lsbValue)
    ecg_fs = ecg.sampleRate
    # rpeaks = movisensobject.getentry('nn_live').values
    # marker = movisensobject.getentry('marker').event
    # artefakte = movisensobject.getentry('artifact').event

    # Plot erstellen

    # ECG
    plt.plot(ecg_signal, label='ECG Signal')

    # # R-Zacken
    # plt.plot(rpeaks[:, 0], ecg_signal[rpeaks[:, 0]],
    #          'r+', label='R-Zacken')
    # # Marker
    # plt.plot(marker, ecg_signal[marker],
    #          'b+', label='Marker')
    # # Artefakte
    # # plt.plot(artefakte, ecg_signal[artefakte],
    # #          'k+', label='Artifacts')

    # Settings
    plt.title('ECG Single Lead movisens ecgMove4 Chest', fontweight="bold")
    plt.xlabel(f'Samples @ {ecg_fs} ')
    plt.ylabel('Magnitude in Bit')

    # Bereiche 10-20 Sekunden
    plt.xlim(10 / (1 / int(ecg_fs)), 20 / (1 / int(ecg_fs)))
    plt.legend()
    plt.show()

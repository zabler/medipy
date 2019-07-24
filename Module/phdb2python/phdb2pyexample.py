import phdb2pyclass
import numpy as np
from matplotlib import pyplot as plt

'''
EXAMPLE

(1) Erstellung eines PHDB Data Objekts mit Signaltyp ECG

(2) ECG Signal, Abtastfrequenz und R-Zacken auslesen

(3) Signal plotten

'''

if __name__ == '__main__':

    # Objekt erstellen, mit Signaltyp ECG
    phdbdata = phdb2pyclass.phdb2pconverter(showtree=False)

    # SignalEntry ECG und ValuesEntry NN_Live auslesen
    ecg_signal = phdbdata.signal
    ecg_fs = phdbdata.fs
    rpeaks = phdbdata.marker
    einheit = phdbdata.units

    # Plot erstellen
    # ECG
    plt.plot(ecg_signal, label='ECG Signal')

    # R-Zacken
    plt.plot(rpeaks, ecg_signal[rpeaks],
             'r+', label='R-Zacken')

    # Settings
    title = phdbdata.filepfad.rsplit('/')[-1] + ' ' + phdbdata.name
    plt.title(f'ECG Single Lead {title}', fontweight="bold")
    plt.xlabel(f'Samples @ {ecg_fs}')
    plt.ylabel(f'Magnitude in {einheit}')

    # Bereiche 10-20 Sekunden
    # plt.xlim(10 / (1 / int(ecg_fs)), 20 / (1 / int(ecg_fs)))
    plt.legend()
    plt.show()

from wda.movesense2python import mo2pclass
from wda.rdetections.dectclass import pantompkins
from matplotlib import pyplot as plt
import numpy as np

'''
MOVESENSE ECG BEATDETECTION

(1) Erstellung eines Movesense Objektes aus ECG Messung

(2) R-Zacken Detektion mit ausgewähltem Algorithmus

(3) Zeitvektor für RR-Intervallplot berechnen

(4) Signale plotten

'''

# Get Data
movesenseobject = mo2pclass.mo2pconverter('ecg')
ecg_signal = movesenseobject.ecg
fs = 128

# Algorithms for Detection
#rpeaks_x, rr = hamilton(ecg_signal, fs).detect()
rpeaks_x, rr = pantompkins(ecg_signal, fs).detect()
# rpeaks_x, rr = ownr(ecg_signal, fs).detect()
# rpeaks_x, rr = skipi(ecg_signal, fs).detect()
# rpeaks_x, rr = peakutils(ecg_signal, fs).detect()

# Calculate Time Vektor for RR
trr = []
trr.append(rr[0])
for k in range(1, len(rr)):
    trr.append(rr[k] + trr[k - 1])

#Plots
plt.figure(1, figsize=(16, 9))
plt.plot(ecg_signal, 'b', label='ECG')
plt.plot(rpeaks_x, ecg_signal[rpeaks_x], 'r+', label='Annotiert')
plt.xlabel(f"Time in Samples @ Sampling frequency {fs} Hz")
plt.ylabel('Amplitude in mV')
plt.title('ECG Signal Single Lead', fontweight="bold")

plt.figure(2, figsize=(16, 9))
plt.plot(trr, rr, 'g')
plt.xlabel(f'Zeit in ms')
plt.ylabel('RR-Intervalle in ms')
plt.title('RR Intervalle', fontweight="bold")

plt.draw()
plt.show()

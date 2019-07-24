from wda.phdb2python import phdb2pyclass
from wda.rdetections.dectclass import pantompkins
from wda.algana import dect_comparer
from matplotlib import pyplot as plt
import numpy as np
import time

'''
MOVESENSE ECG BEATDETECTION

(1) Erstellung eines PHDB Objektes aus ECG Messung und HR

(2) R-Zacken Detektion mit ausgewähltem Algorithmus und Zeiberechnung

(3) Zeitvektor für RR-Intervallplot berechnen

(4) Signale plotten

(5) Algorithmus Bewertung auf Genauigkeit

'''

# Get Data
phdbdata = phdb2pyclass.phdb2pconverter()
ecg_signal = phdbdata.signal
fs = phdbdata.fs
refpeaks = phdbdata.marker

# Algorithm for Detection
start = time.time()
# Check Module/rdetections/detectors for other
#rpeaks_x, rr = hamilton(ecg_signal, fs).detect()
rpeaks_x, rr = pantompkins(ecg_signal, fs).detect()
# rpeaks_x, rr = ownr(ecg_signal, fs).detect()
# rpeaks_x, rr = skipi(ecg_signal, fs).detect()
# rpeaks_x, rr = peakutils(ecg_signal, fs).detect()
end = time.time()
print('________________________________________________________________')
print(f"Algorithm needs {end-start}s for detection\n")

# Calculate Time Vektor for RR
trr = []
trr.append(rr[0])
for k in range(1, len(rr)):
    trr.append(rr[k] + trr[k - 1])

'''Testplot'''
plt.figure(1, figsize=(16, 9))
plt.plot(ecg_signal, 'b', label='ECG')
plt.plot(rpeaks_x, ecg_signal[rpeaks_x], 'r+', label='Annotiert')
plt.xlabel(f'Time in Samples @ Sampling frequency {fs} Hz')
plt.ylabel('Amplitude in mV')
plt.title('ECG Signal Single Lead',
          fontweight="bold")

plt.figure(2, figsize=(16, 9))
plt.plot(trr, rr, 'g')
plt.xlabel(f'Zeit in ms')
plt.ylabel('RR-Intervalle in ms')
plt.title('Rr Intervall time series', fontweight="bold")

# Compare to Reference
dect_comparer.mitcompare(ecg_signal, refpeaks,
                         rpeaks_x, sampleaccuracy=50)
plt.show()

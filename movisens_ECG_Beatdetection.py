from Module.movisens2python import m2pclass
from Module.rdetections.dectclass import*
from Module.algana import dect_comparer
from matplotlib import pyplot as plt
import time
import numpy as np

# Get Data
movisensobject = m2pclass.m2pconverter('ecg')
ecg = movisensobject.getentry('ecg')
ecg_signal = (ecg.signal - int(ecg.baseline)) * float(ecg.lsbValue)
fs = int(ecg.sampleRate)

# Algorithm for Detection
start = time.time()
rpeaks_x, rr = hamilton(ecg_signal, fs).detect()
# rpeaks_x, rr = pantompkins(ecg_signal, fs).detect()
# rpeaks_x, rr = ownr(ecg_signal, fs).detect()
# rpeaks_x, rr = skipi(ecg_signal, fs).detect()
# rpeaks_x, rr = peakutils(ecg_signal, fs).detect()
# Check Module/rdetections/detectors for other
end = time.time()
print('________________________________________________________________')
print(f"Algorithm needs {end-start}s for detection\n")

'''Testplot'''
plt.figure(1, figsize=(16, 9))
plt.plot(ecg_signal, 'b', label='ECG')
plt.plot(rpeaks_x, ecg_signal[rpeaks_x], 'r+', label='Annotiert')
plt.xlabel(f'Time in Samples @ Sampling frequency {fs} Hz')
plt.ylabel('Amplitude in mV')
plt.title('ECG Signal Single Lead',
          fontweight="bold")

# Calculate Time Vektor for RR
trr = []
trr.append(rr[0])
for k in range(1, len(rr)):
    trr.append(rr[k] + trr[k - 1])

plt.figure(2, figsize=(16, 9))
plt.plot(trr, rr, 'g')
plt.xlabel(f'Zeit in ms')
plt.ylabel('RR-Intervalle in ms')
plt.title('Rr Intervall time series', fontweight="bold")

plt.show()

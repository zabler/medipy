from Module.phdb2python import phdb2pyclass
from Module.rdetections.dectclass import*
from Module.algana import dect_comparer
from matplotlib import pyplot as plt
import numpy as np
import time

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

# Get Signalvalues for Rpeaks
rpeaks_y = []
for idx, val in enumerate(ecg_signal):
    if idx in rpeaks_x:
        rpeaks_y.append(val)
rpeaks_y = np.asarray(rpeaks_y)

# Calculate Time Vektor for RR
trr = []
trr.append(rr[0])
for k in range(1, len(rr)):
    trr.append(rr[k] + trr[k - 1])

'''Testplot'''
plt.figure(1, figsize=(16, 9))
plt.plot(ecg_signal, 'b', label='ECG')
plt.plot(rpeaks_x, rpeaks_y, 'r+', label='Annotiert')
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

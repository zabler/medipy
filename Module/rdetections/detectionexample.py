from scipy.misc import electrocardiogram
from dectclass import*
from matplotlib import pyplot as plt
# 'from matplotlib import style
# style.use('fivethirtyeight')'
import numpy as np

'''Testsignal'''
# Hardcoded Testsignal
signal = electrocardiogram()
fs = 360

'''Testdetection'''
# Change Comment for other detectors
rpeaks_x, rr = hamilton(signal, fs).detect()
# rpeaks_x, rr = pantompkins(signal, fs).detect()
# rpeaks_x, rr = ownr(signal, fs).detect()
# rpeaks_x, rr = skipi(signal, fs).detect()
# rpeaks_x, rr = peakutils(signal, fs).detect()

# Get Signalvalues for Rpeaks
rpeaks_y = []
for idx, val in enumerate(signal):
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
plt.plot(signal, 'b', label='ECG')
plt.plot(rpeaks_x, rpeaks_y, 'r+', label='Annotiert')
plt.xlabel(f'Time in Samples @ Sampling frequency {fs} Hz')
plt.ylabel('Amplitude in mV')
plt.title('ECG Signal Single Lead', fontweight="bold")

plt.figure(2, figsize=(16, 9))
plt.plot(trr, rr, 'g')
plt.xlabel(f'Time in ms')
plt.ylabel('RR-Intervalle in ms')
plt.title('Rr Intervall time series', fontweight="bold")

plt.legend()
plt.draw()
plt.show()

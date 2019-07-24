from Module.movesense2python import mo2pclass
from matplotlib import pyplot as plt
import numpy as np

# Get Data
movesenseobject = mo2pclass.mo2pconverter('ecg')
ecg_signal = movesenseobject.ecg
ecg_timestamps = movesenseobject.ecgtime
fs = 128

'''Testplot'''
plt.figure(1, figsize=(16, 9))
plt.plot(ecg_timestamps, ecg_signal, 'b', label='ECG')
plt.xlabel(f"Time in Milis @ Sampling frequency {fs} Hz")
plt.ylabel('Amplitude in µV')
plt.title('ECG Signal Single Lead', fontweight="bold")
plt.draw()
plt.show()

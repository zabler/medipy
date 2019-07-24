from Module.movesense2python import mo2pclass
from Module.rdetections.dectclass import*
from matplotlib import pyplot as plt
import numpy as np
from Module.hrvparameter.hrvclass import hrvparameter
from Module.movesense2python import mo2pclass
from Module.helper import loggertool

# Init Logger
log = loggertool.initlog('movesense_ECG_HRVCalculator')

# Get Data
movesenseobject = mo2pclass.mo2pconverter('ecg')
ecg_signal = movesenseobject.ecg
fs = 128

# Algorithm for Detection
#rpeaks_x, rr = hamilton(ecg_signal, fs).detect()
rpeaks_x, rr = pantompkins(ecg_signal, fs).detect()
# rpeaks_x, rr = ownr(ecg_signal, fs).detect()
# rpeaks_x, rr = skipi(ecg_signal, fs).detect()
# rpeaks_x, rr = peakutils(ecg_signal, fs).detect()

# Berechnung verschiedeneer HRV Parameter
hrv = hrvparameter(rr)
hrv.normalize()
hrv.calcparamter()

# Loggin HRV parameter
for k in hrv.time_domain_features:
    loggertool.log(str(k) + ': ' + str(hrv.time_domain_features[k]), log)
for k in hrv.geometrical_features:
    loggertool.log(str(k) + ': ' + str(hrv.geometrical_features[k]), log)
for k in hrv.frequency_domain_features:
    loggertool.log(str(k) + ': ' + str(hrv.frequency_domain_features[k]), log)
for k in hrv.csi_cvi_features:
    loggertool.log(str(k) + ': ' + str(hrv.csi_cvi_features[k]), log)
for k in hrv.poincare_plot_features:
    loggertool.log(str(k) + ': ' + str(hrv.poincare_plot_features[k]), log)


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
plt.title('ECG Signal Single Lead detected',
          fontweight="bold")

plt.figure(2, figsize=(16, 9))
plt.plot(trr, rr, 'g')
plt.xlabel(f'Zeit in ms')
plt.ylabel('RR-Intervalle in ms')
plt.title('Rr Intervall time series', fontweight="bold")

hrv.plot_distrib()
hrv.plot_poincare()
hrv.plot_psd()
hrv.plot_timeseries()

plt.draw()
plt.show()

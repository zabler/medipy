from wda.hrvparameter.hrvclass import hrvparameter
from wda.helper import loggertool
from wda.movesense2python import mo2pclass
from wda.rdetections.dectclass import pantompkins
from matplotlib import pyplot as plt
import numpy as np


'''
MOVESENSE ECG HRVCALCULATOR

(0) Logger Initiieren

(1) Erstellung eines Movesense Objektes aus ECG Messung

(2) R-Zacken Detektion mit ausgewähltem Algorithmus

(3) Berechnung verschiedener HRV Parameter inkl. Logging

(4) Zeitvektor für RR-Intervallplot berechnen

(5) Signale und Features plotten

'''

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

# Calculate Time Vektor for RR
trr = []
trr.append(rr[0])
for k in range(1, len(rr)):
    trr.append(rr[k] + trr[k - 1])

#Plots
plt.figure(1, figsize=(16, 9))
plt.plot(ecg_signal, 'b', label='ECG')
plt.plot(rpeaks_x, ecg_signal[rpeaks_x], 'r+', label='Annotiert')
plt.xlabel(f'Time in Samples @ Sampling frequency {fs} Hz')
plt.ylabel('Amplitude in mV')
plt.title('ECG Signal Single Lead detected',
          fontweight="bold")

plt.figure(2, figsize=(16, 9))
plt.plot(trr, rr, 'g')
plt.xlabel(f'Zeit in ms')
plt.ylabel('RR-Intervalle in ms')
plt.title('Rr Intervall time series', fontweight="bold")

plt.draw()

hrv.plot_distrib()
hrv.plot_poincare()
hrv.plot_psd()
hrv.plot_timeseries()


plt.show()

from wda.hrvparameter.hrvclass import hrvparameter
from wda.movesense2python import mo2pclass
from wda.helper import loggertool
from matplotlib import pyplot as plt

'''
MOVESENSE HR HRVCALCULATOR

(0) Logger Initiieren

(1) Erstellung eines Movesense Objektes aus HR Messung

(2) Berechnung verschiedener HRV Parameter inkl. Logging

(3) Features plotten

'''

# Init Logger
log = loggertool.initlog('movesense_HR_HRVCalculator')

# Get movesense Data nicht HR sondern Rzacken als Differenzwerte
movesenseobject = mo2pclass.mo2pconverter('heartrate')

# Berechnung verschiedeneer HRV Parameter
hrv = hrvparameter(movesenseobject.rr)
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

hrv.plot_distrib()
hrv.plot_poincare()
hrv.plot_psd()
hrv.plot_timeseries()
plt.show()

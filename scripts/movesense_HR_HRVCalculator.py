from Module.hrvparameter.hrvclass import hrvparameter
from Module.movesense2python import mo2pclass
from matplotlib import pyplot as plt

from Module.helper import loggertool

# Init Logger
log = loggertool.initlog('movesense_HR_HRVCalculator')

# Get movesense Data nicht HR sondern Rzacken als Differenzwerte
movesenseobject = mo2pclass.mo2pconverter('heartrate')
# print(movesenseobject.rr)

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

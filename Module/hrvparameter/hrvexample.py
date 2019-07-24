from scipy.misc import electrocardiogram
from Module.rdetections.dectclass import rdetections
from matplotlib import pyplot as plt
from Module.hrvparameter.hrvclass import hrvparameter
import numpy as np

signal = electrocardiogram()
fs = 360
rpeaks_x, rr = rdetections(signal, fs).hamiltonalg()
hrv = hrvparameter(rr)
hrv.normalize()
hrv.calcparamter()
hrv.plot_distrib()
hrv.plot_poincare()
hrv.plot_psd()
hrv.plot_timeseries()
plt.show()

import hrvanalysis as hrv
import numpy as np

# Only for Example
from scipy.misc import electrocardiogram
from matplotlib import pyplot as plt
from wda.rdetections.dectclass import hamilton


class hrvparameter():
    def __init__(self, rr):
        self.rr = rr
        self.nn = []
        self.time_domain_features = []
        self.geometrical_features = []
        self.frequency_domain_features = []
        self.csi_cvi_features = []
        self.sampen = []

    def normalize(self):
        rr_intervals_without_outliers = hrv.remove_outliers(
            rr_intervals=self.rr, low_rri=300, high_rri=2000)
        interpolated_rr_intervals = hrv.interpolate_nan_values(
            rr_intervals=rr_intervals_without_outliers, interpolation_method="linear")
        nn_intervals_list = hrv.remove_ectopic_beats(
            rr_intervals=interpolated_rr_intervals, method="malik")
        interpolated_nn_intervals = hrv.interpolate_nan_values(
            rr_intervals=nn_intervals_list)
        if np.isnan(interpolated_nn_intervals[-1]):
            interpolated_nn_intervals.pop()
        self.nn = [round(x) for x in interpolated_nn_intervals]

    def calcparamter(self):
        self.time_domain_features = hrv.get_time_domain_features(self.nn)
        self.geometrical_features = hrv.get_geometrical_features(self.nn)
        self.frequency_domain_features = hrv.get_frequency_domain_features(
            self.nn)
        self.csi_cvi_features = hrv.get_csi_cvi_features(self.nn)
        self.poincare_plot_features = hrv.get_poincare_plot_features(self.nn)

    def plot_distrib(self):
        hrv.plot_distrib(self.nn)

    def plot_timeseries(self):
        hrv.plot_timeseries(self.nn)

    def plot_psd(self):
        hrv.plot_psd(self.nn, method="welch")
        hrv.plot_psd(self.nn, method="lomb")

    def plot_poincare(self):
        hrv.plot_poincare(self.nn)

if __name__ == "__main__":
    
    '''
    EXAMPLE

    (1) Erstellung eines PHDB Data Objekts mit Signaltyp ECG

    (2) ECG Signal, Abtastfrequenz und R-Zacken auslesen

    (3) Signal plotten

    '''
    
    signal = electrocardiogram()
    fs = 360
   
    rpeaks_x, rr = hamilton(signal, fs).detect()
    hrvvals = hrvparameter(rr)
    hrvvals.normalize()
    hrvvals.calcparamter()
    hrvvals.plot_distrib()
    hrvvals.plot_poincare()
    hrvvals.plot_psd()
    hrvvals.plot_timeseries()
    plt.show()

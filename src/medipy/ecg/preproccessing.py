'''
preprocessing.py
'''
import scipy as sc
import numpy as np


def hamilton(samples, sample_rate):
    '''
    This method preprocesses a ecg signal for hamilton r peak detection algorithm
    '''
    # Akausaler Bandpass Filter nach [Sedghamiz2014]
    f1 = 5 / sample_rate
    f2 = 15 / sample_rate
    b, a = sc.butter(1, [f1 * 2, f2 * 2], btype='bandpass')
    samples_filtered = sc.filtfilt(b, a, samples)

    # Akusaler 5-Point-Differentiator nach [Pan1985]
    b = [1, 2, 0, -2, -1]
    b = [x * (1 / 8) * self.sample_rate for x in b]
    a = [1]
    samples_diff = sc.filtfilt(b, a, samples_filtered)

    # Betragsbildung (Rectifier)
    samples_rect = abs(samples_diff)

    # Akausaler ungerader 84ms Mittelwertfilter
    b = np.ones(int(0.084 * sample_rate))  # Filterbreite 21 Samples, damit ungerade
    b = b / int(0.084 * sample_rate)
    a = [1]
    samples_ma = sc.filtfilt(b, a, samples_rect)

    return samples_ma

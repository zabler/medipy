'''
r_peak_detection.py
'''

import numpy as np
from scipy import signal
import peakutils as pu


def srinivasan(eeg_raw, eeg_preprocessed, sample_rate):

    # Initiate Return Value
    r_peaks = []

    # Epoch of 1 sec and after 60s check for bullshit
    for i in range(0, len(eeg_preprocessed), sample_rate):
        epoch_1 = eeg_preprocessed[i:i + sample_rate]

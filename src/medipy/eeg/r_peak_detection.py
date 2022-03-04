'''
r_peak_detection.py
'''

import numpy as np


def srinivasan(eeg_preprocessed, sample_rate):
    '''
    This method by srinivasan detects r_peaks of an processed eeg signal
    '''
    # Initiate Return Value
    r_peaks = []

    # Extraction over epoch of 1 sec
    for i in range(0, len(eeg_preprocessed), sample_rate):

        samples_window = eeg_preprocessed[i:i + sample_rate]
        t_val = max(samples_window) * 0.5

        # True R Peak Detection
        for ind, val in enumerate(samples_window):
            if val > t_val:
                r_peaks.append(ind + i)  # i equals the global pos

    return r_peaks

'''
rr_interval_error_correction.py
'''
import math
import numpy as np


def srinivasan(rr_with_nans, sample_rate):
    '''
    This method by srinivasan detects wrong rr intervals and removes them from the list
    '''
    rr_corrected = rr_with_nans
    for i in range(0, len(rr_with_nans), sample_rate):
        rr_window = rr_with_nans[i:i + sample_rate]
        mean_rr = np.nanmean(rr_window)
        for ind, rr_val in enumerate(rr_window):
            if not np.isnan(rr_val):
                if (rr_val < 0.7 * mean_rr) or (rr_val > 1.5 * mean_rr):
                    rr_corrected[ind + i] = np.NaN
    return rr_corrected


#     # Extraction over epoch of 1 sec and after 60s retorspective check for unwanted r-peaks
#     current_rpeaks = []
#     counter = 0
#     for i in range(0, len(eeg_preprocessed), sample_rate):
#         counter += 1
#         samples_window = eeg_preprocessed[i:i + sample_rate]
#         t_val = max(samples_window) * 0.5

#         # True R Peak Detection
#         for ind, val in enumerate(samples_window):
#             if val > t_val:
#                 current_rpeaks.append[ind+i] # i equals the global pos

#         # Retrospective Control
#         if counter == 60:
#            counter = 0
#            mean_rr = np.mean(np.diff(current_rpeaks))
#            for val in current_rpeaks:


#         # r_peak_control()
#         # r_peaks.append()
#         # current_rpeaks = []
    #rr_vals = rr_window[~np.isnan()]
    # rr_val != np.NaN:

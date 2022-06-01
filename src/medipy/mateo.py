'''
mateo.py

Detection of Anomalies and Correction of Ectopic Beats by IPFM Model from Mateo et al. 2003

'''

import numpy as np


def error_detection(r_peaks):
    '''
    Methods A & Appendix
    '''
    def instant_hr(peak_before, peak_current, peak_after):  # Langrange Polynomial
        return 2 * np.abs(np.divide(peak_before - 2 * peak_current + peak_after, (peak_before - peak_current) * (peak_before - peak_after) * (peak_current - peak_after)))

    errors = np.zeros(len(r_peaks))
    threshold_U = 4.3  # * var(r'(5min)), in any case max  at least 0.5 ; nonstat. more than 5min winwdo; assume stationariy signal? then #timeindex +r_peaks? necessary

    for ind, r_peak in enumerate(r_peaks):  # ind+1 = true, becuase we start of initial conidition with the second peak
        # Skip first and last
        if ind not in range(1, len(r_peaks) - 1):
            continue
        elif instant_hr(r_peaks[ind - 1], r_peak, r_peaks[ind + 1]) < threshold_U:
            continue
        elif instant_hr(r_peaks[ind - 1], r_peaks[ind + 1], r_peaks[ind + 2]) < threshold_U:  # Test 1: Remove 0
            # Then FP of 0
            errors[ind] = 1
        elif instant_hr(r_peaks[ind - 1], r_peak, r_peaks[ind + 2]) < threshold_U:  # Test 2: Remove 1
            # Then FP of 1
            errors[ind] = 1
        # elif:  # Test 3: Insert extra between -1 and 0
        #     # Then FN
        # elif:  # Test 4: Insert extra between 0 and 1
        #     # Then FN
        # elif:  # Test 5: Move 0 to Medianpos between -1 and 1
        #     # Then Ectopic
        # elif:  # Test 6: Move 1 to Medianpos between 0 and 2
        #     # Then Ectopic

        # else:  # Consecutive FP, FN or Ectopic
        #     mods = [2, 3, 4, 5]  # number of modifaction to make instant_hr satisfy condition U
        #     # Test 6 conditions abnormal r_peak(0) or r_peaks[ind + 1](1) is errourness
        #     for mod in mods:

    return errors  # dict


def error_correction():
    '''
    Methods B
    '''
    n_peaks = 404
    #!Delete extras/fp!
    # misseds/fn and ectopics?
    return n_peaks  # dict

'''
rr_interval_calculation.py
'''
import math


def calculator(r_peaks):
    '''
    This method calculates all rr intervals of a given r_peak list
    '''
    rr_intervals = []
    for k in range(1, len(r_peaks)):
        rr_intervals.append(math.ceil(r_peaks[k] - r_peaks[k - 1]))
    return rr_intervals

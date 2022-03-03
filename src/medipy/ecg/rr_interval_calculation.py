'''
rr_interval_calculation.py
'''
import math


def calculator(r_peaks, sample_rate=None):
    '''
    This method calculates all rr intervals of a given r_peak list
    '''
    period_ms = int((1 / sample_rate) * 1000)
    rr_intervals = []
    for k in range(1, len(r_peaks)):
        rr_intervals.append(math.ceil(r_peaks[k] - r_peaks[k - 1]))
    # an Grid anpassen
    if sample_rate is not None:
        return [rr_interval * period_ms for rr_interval in rr_intervals]
    else:
        return rr_intervals

'''
hrv_feature_calcualtion.py
'''

import numpy as np
import nolds
from astropy.timeseries import LombScargle


def hrv_features_time(rr_intervals):
    '''
    This method calculates all time features
    '''

    # For Calculation
    drr = np.diff(rr_intervals)
    heart_rate = np.divide(60000, rr_intervals)

    # Statistical Features
    sdnn = np.std(rr_intervals, ddof=1)
    sdsd = np.std(drr, ddof=1)
    rmssd = np.sqrt(np.mean(drr ** 2))
    nn50 = sum(np.abs(drr) > 50)
    pnn50 = (nn50 / len(rr_intervals)) * 100

    # Heartrate Features
    hr_max_min = max(heart_rate) - min(heart_rate)

    time_features = {
        'SDNN': sdnn,
        'SDSD': sdsd,
        'RMSSD': rmssd,
        'NN50': nn50,
        'pNN50': pnn50,
        'HR_MAX_MIN': hr_max_min,
    }
    return time_features


def hrv_features_frequency(rr_intervals):
    '''
    This method calculates all frequency features
    '''
    # For Calculation
    rr_timestamps_cumsum = np.cumsum(rr_intervals) / 1000  # in sec damit Hz
    rr_timestamps = rr_timestamps_cumsum - rr_timestamps_cumsum[0]

    # LombScargle by Astropy
    freq, psd = LombScargle(rr_timestamps, rr_intervals, normalization='psd').autopower(minimum_frequency=0.040, maximum_frequency=0.400, samples_per_peak=5)

    # LF Band
    lf_indexes = np.logical_and(freq >= 0.04, freq < 0.15)

    # HF Band
    hf_indexes = np.logical_and(freq >= 0.15, freq < 0.40)

    # Spectal Features
    lf_power = np.trapz(y=psd[lf_indexes], x=freq[lf_indexes])
    hf_power = np.trapz(y=psd[hf_indexes], x=freq[hf_indexes])
    lf_peak = freq[lf_indexes].max()
    hf_peak = freq[hf_indexes].max()
    lf_hf_ratio = lf_power / hf_power
    lf_nu = (lf_power / (lf_power + hf_power)) * 100
    hf_nu = (hf_power / (lf_power + hf_power)) * 100

    freqency_features = {
        'LF_POWER': lf_power,
        'HF_POWER': hf_power,
        'LF_HF_RATIO': lf_hf_ratio,
        'LF_NU': lf_nu,
        'HF_NU': hf_nu,
        'LF_PEAK': lf_peak,
        'HF_PEAK': hf_peak
    }

    return freqency_features


def hrv_features_nonlinear(rr_intervals):
    '''
    This method calculates all nonlinear features
    '''
    # For Calculation
    drr = np.diff(rr_intervals)
    short = range(4, 16 + 1)
    long = range(17, 64 + 1)

    # Poincare Features
    sd1 = np.sqrt((np.std(drr, ddof=1) ** 2) * 0.5)
    sd2 = np.sqrt(2 * np.std(rr_intervals, ddof=1) ** 2 - 0.5 * np.std(drr, ddof=1) ** 2)
    t_feature = 4 * sd1
    l_feature = 4 * sd2
    csi = sd2 / sd1
    modified_csi = (l_feature ** 2) / t_feature
    cvi = np.log10(l_feature * t_feature)

    # DFA Features
    if len(rr_intervals) > 2:
        df_alpha_1 = nolds.dfa(rr_intervals, short, debug_data=False, overlap=False)
        df_alpha_2 = nolds.dfa(rr_intervals, long, debug_data=False, overlap=False)
    else:
        df_alpha_1 = np.nan
        df_alpha_2 = np.nan

    nonlinear_features = {
        'SD1': sd1,
        'SD2': sd2,
        'CSI': csi,
        'MODIFIED_CSI': modified_csi,
        'CVI': cvi,
        'DF_ALPHA_1': df_alpha_1,
        'DF_ALPHA_2': df_alpha_2
    }

    return nonlinear_features

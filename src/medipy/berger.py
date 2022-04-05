'''
berger.py

Berger algorithm implemented in Python by Nicolas Zabler, 2022.

This file contains an algorithm, developed by [1] in 1986 to resample and interpolate an non-equidistant list of r-peaks, extracted from an ECG into a RR intervall time series.
According to [2], the algorithm is goo.. this algorithm is better than other

Input
rr: List of RR intervals in ms
resample_fs: Resample Frequency

Output
Vector with equidistant, resampled and interpolated RR interval values

Recommendations
Use a resample frequency which is a multiple of the origin sample frequency, so the resulting interpolated RR array is easy to contenate

References
[1] R. D. Berger, S. Akselrod, D. Gordon and R. J. Cohen, "An Efficient Algorithm for Spectral Analysis of Heart Rate Variability," in IEEE Transactions on Biomedical Engineering, vol. BME-33, no. 9, pp. 900-904, Sept. 1986, doi: 10.1109/TBME.1986.325789.
[2] J. Vila, S. Barro, J. Presedo, R. Ruiz and F. Palacios, "Analysis of heart rate variability with evenly spaced time values," 1992 14th Annual International Conference of the IEEE Engineering in Medicine and Biology Society, 1992, pp. 575-576, doi: 10.1109/IEMBS.1992.5761116.
[3] Singh D, Vinod K, Saxena SC. Sampling frequency of the RR interval time series for spectral analysis of heart rate variability. J Med Eng Technol. 2004 Nov-Dec;28(6):263-72. doi: 10.1080/03091900410001662350. PMID: 15513744.
'''

import numpy as np


def resampling(rr, resample_fs=4):
    # Init Resampling
    period_ms = int(1000 / resample_fs)
    r_peaks = np.cumsum(rr)
    r_peaks = np.insert(r_peaks, 0, 0)
    grid = np.arange(0, r_peaks[-1], period_ms)  # Not sure if correct
    rr_interpolated = []  # change to vector with know length
    rr.append(rr[-1])
    rr_current = 0
    # rr_interpolated = np.empty(int(np.divide(np.sum(rr), period_ms)))  # Vector length = Number of Vals = TotalTime(sum(rr)) * resample_fs  = TotalTime(sum(rr)) / period_ms

    # Iterating over grid, starting with first R-peak, ends with or before oder last R-Peak, depending on origin fs
    for val in grid:  # grid[1:]
        r_peaks_local = np.where(np.logical_and(r_peaks >= val - period_ms, r_peaks < val + period_ms))[0]
        if r_peaks_local.size == 0:
            rri_int = rr[rr_current]  # CHANGE? # Look for closest r-peaks with quick algo?
        else:
            rr_vals_local = rr[int(min(r_peaks_local)):int(max(r_peaks_local)) + 2]  # if more than one r-peak, mid RR-values weighted 1 by formular automatically
            rri_int = int(np.divide(2 * period_ms, (np.divide(r_peaks[r_peaks_local[0]] - (val - period_ms), rr_vals_local[0]) +
                                                    np.divide(val + period_ms - r_peaks[r_peaks_local[-1]], rr_vals_local[-1]))))
            rr_current += 1  # not correct ?! at least doest not work
        rr_interpolated.append(rri_int)

    return rr_interpolated

# Mögliche Test
# Is fs vielfaches von rfs?
# Wenn nein quit 7 oder Warnung # Modulo Operator
# Evtl. als function implementieren
# rri_int = _calc_rri_int(val, period_ms, list(r_peaks_local), rr_vals_local)
# def _calc_rri_int(t_0, delta_t, local_r_peaks, local_rr_vals):

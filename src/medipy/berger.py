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
First R-peak of signal equals the first value of returend interpolated RR List

References
[1] R. D. Berger, S. Akselrod, D. Gordon and R. J. Cohen, "An Efficient Algorithm for Spectral Analysis of Heart Rate Variability," in IEEE Transactions on Biomedical Engineering, vol. BME-33, no. 9, pp. 900-904, Sept. 1986, doi: 10.1109/TBME.1986.325789.
[2] J. Vila, S. Barro, J. Presedo, R. Ruiz and F. Palacios, "Analysis of heart rate variability with evenly spaced time values," 1992 14th Annual International Conference of the IEEE Engineering in Medicine and Biology Society, 1992, pp. 575-576, doi: 10.1109/IEMBS.1992.5761116.
[3] Singh D, Vinod K, Saxena SC. Sampling frequency of the RR interval time series for spectral analysis of heart rate variability. J Med Eng Technol. 2004 Nov-Dec;28(6):263-72. doi: 10.1080/03091900410001662350. PMID: 15513744.
'''

import numpy as np


def resampling(rr, resample_fs=4):
    # Init Resampling
    rr = list(rr)
    period_ms = int(1000 / resample_fs)
    r_peaks = np.cumsum(rr)
    r_peaks = np.insert(r_peaks, 0, 0)
    grid = np.arange(0, r_peaks[-1], period_ms)
    rr_interpolated = np.empty(len(grid))
    rr.append(rr[-1])
    rr_current = rr[0]

    # First value on t_0 to RR_0
    rr_interpolated[0] = rr[0]

    # Iterating over grid, starting with t_1 after first R-peak(t_0), ends with or before oder last R-Peak, depending on origin fs
    for ind, val in enumerate(grid[1:]):  #
        r_peaks_local = np.where(np.logical_and(r_peaks > (val - period_ms), r_peaks <= (val + period_ms)))[0]
        if r_peaks_local.size == 0:
            rri_int = rr_current
        else:
            rr_vals_local = rr[int(min(r_peaks_local)) - 1:int(max(r_peaks_local)) + 1]  # if more than one r-peak, mid RR-values weighted 1 by formular automatically
            rri_int = int(np.divide(2 * period_ms, (np.divide(r_peaks[r_peaks_local[0]] - (val - period_ms), rr_vals_local[0]) +
                                                    np.divide(val + period_ms - r_peaks[r_peaks_local[-1]], rr_vals_local[-1]))))
            rr_current = rr_vals_local[-1]
        rr_interpolated[ind + 1] = rri_int
    return rr_interpolated

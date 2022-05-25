'''
lazaro.py

Lazaro algorithm implemented in Python by Nicolas Zabler, 2022.

This file contains an algorithm, developed by [1] in 2014 extract a respiratory signal with the means of slopes in an QRS complex.
That pattern was first found by [2].

Input
N-peaks: list of lists: [cleaned R-peak, corresponding Q_peak, corresponding S_peak]
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

import math
import numpy as np


def slope_extraction(ecg, n_peaks, q_n_peaks, s_n_peaks, resample_frequency):

    # TO-DO
    # Generalize with resample frequency at LSQ: np.linalg.lstsq(A, y, rcond=None)[0] to always use 8ms windows accoringt to lazaro or romero 15ms
    # Not neceassary for current approach because of 250Hz sample rate which makes 4ms windows
    # Implement 2nd search
    # make least sqaure not mean

    # Init
    up_slope_vector = np.empty(len(n_peaks))
    up_slope_vector[:] = np.NaN
    down_slope_vector = np.empty(len(n_peaks))
    down_slope_vector[:] = np.NaN
    peak_angle_vector = np.empty(len(n_peaks))
    peak_angle_vector[:] = np.NaN
    resample_frequency = int(resample_frequency)

    # Calc
    for ind, n_peak in enumerate(n_peaks):
        # Check if q and s are availabel
        if math.isnan(q_n_peaks[ind]):
            continue
        else:
            q_n = int(q_n_peaks[ind])
        if math.isnan(s_n_peaks[ind]):
            continue
        else:
            s_n = int(s_n_peaks[ind])

        # Steigung des Least Square Fit +/- 4ms um Slope_Point (=+/- Sample bei 250Hz) (1 Sample entspricht 4ms, theoretisch reicht np. diff max aus; jetzt aber LS um 3 Werte;
        linear_matrix = np.vstack([[0, 1, 2], [1, 1, 1]]).T

        up_slope_point = np.argmax(np.diff(ecg[q_n:n_peak]))
        up_slope, _ = np.linalg.lstsq(linear_matrix, ecg[q_n + up_slope_point - 1:q_n + up_slope_point + 2], rcond=None)[0]
        up_slope = np.max(np.diff(ecg[q_n:n_peak]))
        up_slope_vector[ind] = up_slope

        down_slope_point = np.argmax(np.diff(ecg[n_peak:s_n]))
        down_slope, _ = np.linalg.lstsq(linear_matrix, ecg[n_peak + down_slope_point - 1:n_peak + down_slope_point + 2], rcond=None)[0]
        down_slope = np.max(np.diff(ecg[n_peak:s_n]))
        down_slope_vector[ind] = down_slope

        peak_angle = np.arctan(np.divide((up_slope - down_slope), 0.4 * (6.25 + (up_slope * down_slope))))
        peak_angle_vector[ind] = peak_angle

    slopes = {'up_slope': up_slope_vector, 'down_slope': down_slope_vector, 'peak_angle': peak_angle_vector}
    return slopes

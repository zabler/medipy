'''
bailon.py

Interpolation by IPFM Model (Interpolate N-Peaks and Differentiate them) from Hernando et al 2016 and Bailon et al 2011

'''
from scipy.interpolate import CubicSpline
import numpy as np
from scipy.signal import kaiserord, lfilter, firwin, freqz


def tvipfm_interpolation(n_peaks_digits, sample_rate=250):
    '''
    Time varying ipfm modualtion ...

    n_peaks list starting from first peak, equals t_0, to last sample of the signal, non peaks are nans

    # the interpolation has the same frequency as the ecg signal itself

    numpy array n_peaks_y with Nans at the specific positiions, so an interp is not needed
    '''

    # N-Peak Pairs
    n_peaks_x = n_peaks_digits  # equals t, digital
    n_peaks_y = np.arange(len(n_peaks_digits))  # equals k(t), digital

    # Build Interpolation Model using (t, k(t)) Pairs using CubicScpline; in Bailon et al. 5 5th order spline is used
    interp_model = CubicSpline(n_peaks_x, n_peaks_y, bc_type='natural', extrapolate=None)

    # Interpolation
    n_peaks_x_new = np.arange(n_peaks_digits[0], n_peaks_digits[-1])  # equals t, continously (same fs as the original signal)
    n_peaks_y_new = interp_model(n_peaks_x_new)  # interpolated kappa, continously (same fs as the original signal)

    # Estimation d_HR
    d_HR = np.diff(n_peaks_y_new)  # equals derived kappa

    # Estimation d_HRM
    # Gruppenlaufzeit konstant; Verzögerung? Subtraction ?; besonders schauen wie stark die Dämpfung für für 0,04-0,4 Hz ist

    # Filter design FILTFILT, FILTER RESPONSE ETC..
    nyq_rate = sample_rate / 2.0  # Filter Settings
    width = 5.0 / nyq_rate  # The desired width of the transition from pass to stop, relative to the Nyquist rate
    ripple_db = 60.0  # The desired attenuation in the stop band, in dB.
    N, beta = kaiserord(ripple_db, width)  # Compute the order and Kaiser parameter for the FIR filter.
    cutoff_hz = 0.03  # The cutoff frequency of the filter.
    taps = firwin(N, cutoff_hz / nyq_rate, window=('kaiser', beta))  # Use firwin with a Kaiser window to create a lowpass FIR filter.
    d_HRM = lfilter(taps, 1.0, d_HR)  # Schätzung von d_HRM durch Anwendung 0,03 Hz LP FIR Filter auf d_HR

    # HRV signal d_HRV UM ZWEI KÜRZER ALS INPUT
    d_HRV = d_HR - d_HRM

    return d_HRV

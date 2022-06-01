'''
bailon.py

Interpolation and Resampling by IPFM Model (Interpolate N-Peaks and Differentiate them) from Hernando et al 2016 and Bailon et al 2011

'''
from scipy.interpolate import CubicSpline
import numpy as np


def tvipfm_interpolation(n_peaks, interp_freq=4):
    '''
    Time varying ipfm modualtion ...
    '''

    # input values, n_peaks_x, n_peaks_y with Nans at the specific positiions, so an interp is not needed

    # N-Peak Pairs
    n_peaks_x = np.arange(len(n_peaks))
    n_peaks_y = n_peaks

    # Interpolation Model # (add time varing threshold add error correction)
    interp_model = CubicSpline(n_peaks_x, n_peaks_y, bc_type='natural', extrapolate=None)

    # New Support Points alle 1/fs ms?? How missing vals, how do i know that a val is missing, startend mit dem ersten N-Peak und endend mit dem letzten Peak/ oder danach
    n_peaks_x_new = np.linspace(0, n_peaks_x[-1], int((1 / interp_freq) * 1000))  # (t_N - t_0) /(1/fs)

    # Interpolated Kappa
    n_peaks_y_new = interp_model(n_peaks_x_new)

    # Estimation d_HR
    # Ableitung bilden / Kappa'

    # Estimation d_HRM
    # LP bilden von d_HR

    # d_HRV
    # d_HRV = d_HR - d_HRM

    return n_peaks_y_new

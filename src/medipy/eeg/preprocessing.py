'''
preprocessing.py
'''
from scipy.fftpack import cc_diff
import scipy.signal as sc
import numpy as np
import pywt


def srinivasan(samples, sample_rate):
    '''
    This method by srinivasan2015 prepares an eeg signal for r-peak detection using (1) Wavelet Decomposition and (2) TK Energy Operator
    '''
    # ana_len = samplfor every s (250 vals)
    samples_window = [2, 3, 4, 2, 2, 3]
    # normalize signal ?!?!?
    samples_window_norm = 0.5 * samples_window
    # the wavelet 2^3 equals 47-63 bandpass because here minimum interference from EEG and ECG
    coeffs = pywt.wavedec(samples_window_norm, 'coif1', level=3)
    cA3, cD3, cD2, cD1 = coeffs
    # What to do with the coeffs? Wie erhalte ich daraus ein Signal zurück?Faltung mit downsamling mehrfach?
    # hochpässe filtern und tiefpässe und die Deetails signlas vom Orignal subtrahieren ?
    # to reduce noise signal is filtered with an teager-kaiser energy operator
    samples_tk_fitlered = cA3

    return samples_tk_fitlered

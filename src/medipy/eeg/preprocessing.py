'''
preprocessing.py
'''
from distutils.log import log
import scipy.signal as sc
import numpy as np
import pywt
from matplotlib import pyplot as plt


def srinivasan(samples, sample_rate):
    '''
    This method by srinivasan2015 prepares an eeg signal for r-peak detection using (1) Wavelet Decomposition and (2) TK Energy Operator
    '''
    # Initiate Return Value
    samples_filtered = []

    # Epoch of 1 sec
    for i in range(0, len(samples), sample_rate):
        samples_window = samples[i:i + sample_rate]

        # Normalize? by what ? mean? energy?
        # Normalisierung BETRAG hilfreich?
        samples_window_norm = np.divide(samples_window - np.min(samples_window), (np.max(samples_window) - np.min(samples_window)))
        # plt.figure()
        # plt.plot(samples_window_norm)

        # DWT Wavelet Decomposition Level 3 with Coiflet 1
        coeffs = pywt.wavedec(samples_window_norm, 'coif1', level=3)
        cA3, cD3, cD2, cD1 = coeffs
        # plt.figure()
        # ax1 = plt.subplot(141)
        # ax1.plot(cA3)
        # ax2 = plt.subplot(142)
        # ax2.plot(cD3)
        # ax3 = plt.subplot(143)
        # ax3.plot(cD2)
        # ax4 = plt.subplot(144)
        # ax4.plot(cD1)

        # USING cA3 vs Thresholding & Wavelet Reconstruction

        N = len(samples_window_norm)
        # Soft Thressholding / True Denoising
        N = len(samples_window_norm)
        t_cD1 = np.divide(np.median(abs(cD1)) * np.sqrt(2 * np.log(len(cD1))), 0.6745)
        cD1_thresh = pywt.threshold(cD1, t_cD1, mode='soft')
        t_cD2 = np.divide(np.median(abs(cD2)) * np.sqrt(2 * np.log(len(cD2))), 0.6745)
        cD2_thresh = pywt.threshold(cD2, t_cD2, mode='soft')
        t_cD3 = np.divide(np.median(abs(cD3)) * np.sqrt(2 * np.log(len(cD3))), 0.6745)
        cD3_thresh = pywt.threshold(cD3, t_cD3, mode='soft')

        #cA3_tresh = pywt.threshold(cA3, tresh_val, mode=mode_val)
        # plt.figure()
        # ax1 = plt.subplot(141)
        # ax1.plot(cA3)
        # ax2 = plt.subplot(142)
        # ax2.plot(cD3_thresh)
        # ax3 = plt.subplot(143)
        # ax3.plot(cD2_thresh)
        # ax4 = plt.subplot(144)
        # ax4.plot(cD1_thresh)

        # Wavelet Reconstruciton
        coeffs_denoised = cA3, cD3_thresh, cD2_thresh, cD1_thresh
        samples_denoised = pywt.waverec(coeffs_denoised, 'coif1')
        # plt.figure()
        # plt.plot(samples_denoised)

        # TK Energy Operator
        samples_filtered_tk = []
        for index, val in enumerate(samples_denoised):
            if (index == 0) or (index == len(samples_denoised) - 1):
                continue
            else:
                samples_filtered_tk.append((val * val) - (samples_denoised[index - 1] * samples_denoised[index + 1]))

        samples_filtered_tk.insert(0, samples_filtered_tk[0])
        samples_filtered_tk.append(samples_filtered_tk[-1])
        # samples_filtered_tk pos 0 und -1 setzen
        # plt.figure()
        # plt.plot(samples_filtered_tk)

        samples_filtered += samples_filtered_tk  # samples_filtered += list(cA3)

        # Extract R peaks here ?

    return samples_filtered

'''
srinivasan.py
'''
import scipy.signal as sc
import math
import numpy as np
import pywt
from matplotlib import pyplot as plt
#from .extlib.modwtpy import modwt

# Class srinivasan?


def preprocessing(samples, sample_rate, ecg, plot=False):
    '''
    This method by srinivasan2015 prepares an eeg signal for r-peak detection using (1) Wavelet Decomposition and (2) TK Energy Operator
    '''
    # Initiate Return Value
    samples_filtered = []
    samples_normed = []

    # Epoch of 1 sec
    for i in range(0, len(samples), sample_rate):
        ecg_window = ecg[i:i + sample_rate]
        samples_window = samples[i:i + sample_rate]
        # Normalization (by what mean? energy? standardisation?)
        samples_window_norm = np.divide(samples_window - np.min(samples_window), (np.max(samples_window) - np.min(samples_window)))
        samples_normed += list(samples_window_norm)
        # samples_window_stand = np.divide(samples_window - np.mean(samples_window), np.std(samples_window))

        if plot:

            fig_pipe = plt.figure()
            fig_wave = plt.figure()

            # # MODWT Wavelet Decomposition Level 3 with Coiflet 1 (32-64) 47-63?? #brauch sehr lange
            # wt = modwt(samples_window_norm, 'coif1', 3)
            # cA3 = wt[3]
            # cD3 = wt[2]  # Decompoistion Output corresponding to 47-63 Hz
            # cD2 = wt[1]
            # cD1 = wt[0]

            # DWT Wavelet Decomposition Level 3 with Coiflet 1 using PyWavelets
            coeffs = pywt.wavedec(samples_window_norm, 'coif1', level=3)
            cA3 = coeffs[0]
            cD3 = coeffs[1]  # Decompoistion Output corresponding to 47-63 Hz
            cD2 = coeffs[2]
            cD1 = coeffs[3]

            # Thressholding / True Denoising
            thresh_mode = 'hard'
            # t_cD1 = np.divide(np.median(abs(cD1)) * np.sqrt(2 * np.log(len(cD1))), 0.6745)
            # cD1_thresh = pywt.threshold(cD1, t_cD1, mode=thresh_mode)
            cD1_thresh = np.zeros(len(cD1))
            # t_cD2 = np.divide(np.median(abs(cD2)) * np.sqrt(2 * np.log(len(cD2))), 0.6745)
            # cD2_thresh = pywt.threshold(cD2, t_cD2, mode=thresh_mode)
            cD2_thresh = np.zeros(len(cD2))
            # t_cD3 = np.divide(np.median(abs(cD3)) * np.sqrt(2 * np.log(len(cD3))), 0.6745)
            # cD3_thresh = pywt.threshold(cD3, t_cD3, mode=thresh_mode)
            # cD3_thresh = np.zeros(len(cD3))
            cD3_thresh = cD3
            # t_cA3 = np.divide(np.median(abs(cA3)) * np.sqrt(2 * np.log(len(cA3))), 0.6745)
            # cA3_thresh = pywt.threshold(cA3, t_cA3, mode=thresh_mode)
            cA3_thresh = np.zeros(len(cA3))

            # # Akausaler Bandpass Filter nach [Sedghamiz2014] #evtl. nur Highpass damit die Baseline nicht so drifftet
            # f1 = 50 / sample_rate  # 5
            # f2 = 60 / sample_rate  # 15
            # b, a = sc.butter(1, [f1 * 2, f2 * 2], btype='bandpass')
            # # #f_limit = 1 / eeg_fs
            # # #b, a = sc.butter(1, f_limit, btype='high')
            # cD3_thresh = sc.filtfilt(b, a, cD3_thresh)

            # Wavelet Reconstruciton
            coeffs_denoised = cA3_thresh, cD3_thresh, cD2_thresh, cD1_thresh
            samples_denoised = pywt.waverec(coeffs_denoised, 'coif1', mode='zero')

            # TK Energy Operator
            samples_filtered_tk = []
            for index, val in enumerate(samples_denoised[1:-1]):
                samples_filtered_tk.append((val * val) - (samples_denoised[index - 1] * samples_denoised[index + 1]))
            samples_filtered_tk.insert(0, samples_filtered_tk[0])
            samples_filtered_tk.append(samples_filtered_tk[-1])

            # Complete List
            samples_filtered += samples_filtered_tk  # samples_filtered += list(cA3)

            # fig_pipe = plt.figure()
            ax10 = fig_pipe.add_subplot(611)
            ax10.plot(ecg_window)
            ax11 = fig_pipe.add_subplot(612)
            ax11.plot(samples_window)
            ax12 = fig_pipe.add_subplot(613)
            ax12.plot(samples_window_norm)
            ax13 = fig_pipe.add_subplot(614)
            ax13.plot(samples_denoised)
            ax14 = fig_pipe.add_subplot(615)
            ax14.plot(samples_filtered_tk)
            ax15 = fig_pipe.add_subplot(616)
            ax15.plot(samples_filtered)

            # #fig_wave = plt.figure()
            ax21 = fig_wave.add_subplot(421)
            ax21.plot(cD1)
            ax22 = fig_wave.add_subplot(423)
            ax22.plot(cD2)
            ax23 = fig_wave.add_subplot(425)
            ax23.plot(cD3)
            ax24 = fig_wave.add_subplot(427)
            ax24.plot(cA3)
            ax25 = fig_wave.add_subplot(422)
            ax25.plot(cD1_thresh)
            ax26 = fig_wave.add_subplot(424)
            ax26.plot(cD2_thresh)
            ax27 = fig_wave.add_subplot(426)
            ax27.plot(cD3_thresh)
            ax28 = fig_wave.add_subplot(428)
            ax28.plot(cA3_thresh)

            plt.ion()
            plt.show()  # block=False)
            # plt.pause(5)
            fig_pipe.clear()
            fig_wave.clear()

            # pass
            plt.close('all')

    # DWT Wavelet Decomposition Level 3 with Coiflet 1 using PyWavelets
    coeffs = pywt.wavedec(samples_normed, 'coif1', level=3)
    cA3 = coeffs[0]
    cD3 = coeffs[1]  # Decompoistion Output corresponding to 47-63 Hz
    cD2 = coeffs[2]
    cD1 = coeffs[3]

    # Thressholding / True Denoising
    thresh_mode = 'hard'
    # t_cD1 = np.divide(np.median(abs(cD1)) * np.sqrt(2 * np.log(len(cD1))), 0.6745)
    # cD1_thresh = pywt.threshold(cD1, t_cD1, mode=thresh_mode)
    cD1_thresh = np.zeros(len(cD1))
    # t_cD2 = np.divide(np.median(abs(cD2)) * np.sqrt(2 * np.log(len(cD2))), 0.6745)
    # cD2_thresh = pywt.threshold(cD2, t_cD2, mode=thresh_mode)
    cD2_thresh = np.zeros(len(cD2))
    # t_cD3 = np.divide(np.median(abs(cD3)) * np.sqrt(2 * np.log(len(cD3))), 0.6745)
    # cD3_thresh = pywt.threshold(cD3, t_cD3, mode=thresh_mode)
    # cD3_thresh = np.zeros(len(cD3))
    cD3_thresh = cD3
    # t_cA3 = np.divide(np.median(abs(cA3)) * np.sqrt(2 * np.log(len(cA3))), 0.6745)
    # cA3_thresh = pywt.threshold(cA3, t_cA3, mode=thresh_mode)
    cA3_thresh = np.zeros(len(cA3))

    # Akausaler Bandpass Filter nach [Sedghamiz2014] #evtl. nur Highpass damit die Baseline nicht so drifftet
    # f1 = 50 / sample_rate  # 5
    # f2 = 60 / sample_rate  # 15
    # b, a = sc.butter(1, [f1 * 2, f2 * 2], btype='bandpass')
    # #f_limit = 1 / eeg_fs
    # #b, a = sc.butter(1, f_limit, btype='high')
    # cD3_thresh = sc.filtfilt(b, a, cD3)

    # Wavelet Reconstruciton
    coeffs_denoised = cA3_thresh, cD3_thresh, cD2_thresh, cD1_thresh
    samples_denoised = pywt.waverec(coeffs_denoised, 'coif1', mode='zero')

    # TK Energy Operator
    samples_filtered_tk = []
    for index, val in enumerate(samples_denoised[1:-1]):
        samples_filtered_tk.append((val * val) - (samples_denoised[index - 1] * samples_denoised[index + 1]))
    samples_filtered_tk.insert(0, samples_filtered_tk[0])
    samples_filtered_tk.append(samples_filtered_tk[-1])

    return samples_filtered_tk


def r_peak_detection(eeg_preprocessed, sample_rate):
    '''
    This method by srinivasan detects r_peaks of an processed eeg signal
    '''
    # Initiate Return Value
    local_peaks = []

    # Extraction over epoch of 1 sec
    for i in range(0, len(eeg_preprocessed), sample_rate):

        samples_window = eeg_preprocessed[i: i + sample_rate]
        t_val = max(samples_window) * 0.9  # 0.5 Original Val

        # Local Peaks
        for ind, val in enumerate(samples_window):
            if val > t_val:
                local_peaks.append(ind + i)

    return local_peaks


def r_peak_correction(eeg_preprocessed, r_peaks, sample_rate):
    mask = np.zeros(len(r_peaks))
    for ind in range(0, len(r_peaks) - 1):
        if (ind == 0) or (ind == len(r_peaks) - 1):
            continue
        if (abs(r_peaks[ind - 1] - r_peaks[ind]) >= 125) and (abs(r_peaks[ind + 1] - r_peaks[ind]) >= 125):  # CHANGE to automatic
            mask[ind] = 1  # Andere Werte sind fern von peak
        else:  # Andere Werte sind nah an peak
            if (eeg_preprocessed[r_peaks[ind - 1]] > eeg_preprocessed[r_peaks[ind]]) or (eeg_preprocessed[r_peaks[ind + 1]] > eeg_preprocessed[r_peaks[ind]]):
                mask[ind] = 0  # Andere Werte sind größer als peak
            else:
                mask[ind] = 1  # Andere Werte sind kleiner als peak
    if mask[1] == 0:  # Ausnahme Erster Wert abhängig von Zweitem Wert, sonst Null lassen
        mask[0] = 1
    if mask[-2] == 0:  # Ausnahme Letzter Wert abhängig von vorletztem Wert, sonst Null lassen
        mask[-1] = 1
    r_peaks = np.array(r_peaks)
    mask_bool = mask.astype(bool)
    r_peaks_filtered = r_peaks[mask_bool]  # Update local peaks list

    return r_peaks_filtered


def rr_interval_error_correction(rr_with_nans, sample_rate):
    '''
    This method by srinivasan detects wrong rr intervals and removes them from the list
    '''
    rr_corrected = rr_with_nans
    for i in range(0, len(rr_with_nans), sample_rate * 60):  # damit auf eine Minute gewindowed ohne Overlapp
        rr_window = rr_with_nans[i:i + sample_rate * 60]
        mean_rr = np.nanmean(rr_window)
        for ind, rr_val in enumerate(rr_window):
            if not np.isnan(rr_val):
                if (rr_val < 0.7 * mean_rr) or (rr_val > 1.5 * mean_rr):
                    rr_corrected[ind + i] = np.NaN
    return rr_corrected


# def preprocessing(samples, sample_rate, ecg, plot=False):
#     '''
#     This method by srinivasan2015 prepares an eeg signal for r-peak detection using (1) Wavelet Decomposition and (2) TK Energy Operator
#     '''

#     # Normalization (by what mean? energy? standardisation?)
#     #samples_norm = np.divide(samples- np.min(samples), (np.max(samples) - np.min(samples)))
#     samples_stand = np.divide(samples - np.mean(samples), np.std(samples))

#     # DWT Wavelet Decomposition Level 3 with Coiflet 1 using PyWavelets
#     coeffs = pywt.wavedec(samples_stand, 'coif1', level=3)
#     cA3 = coeffs[0]
#     cD3 = coeffs[1]  # Decompoistion Output corresponding to 47-63 Hz
#     cD2 = coeffs[2]
#     cD1 = coeffs[3]

#     # # Akausaler Bandpass Filter nach [Sedghamiz2014] #evtl. nur Highpass damit die Baseline nicht so drifftet
#     # f1 = 1 / eeg_fs  # 5
#     # f2 = 70 / eeg_fs  # 15
#     # b, a = sc.butter(1, [f1 * 2, f2 * 2], btype='bandpass')
#     # #f_limit = 1 / eeg_fs
#     # #b, a = sc.butter(1, f_limit, btype='high')
#     # eeg_preprocessed = sc.filtfilt(b, a, eeg_raw)

#     # Thressholding / True Denoising
#     thresh_mode = 'hard'
#     t_cD1 = np.divide(np.median(abs(cD1)) * np.sqrt(2 * np.log(len(cD1))), 0.6745)
#     #cD1_thresh = pywt.threshold(cD1, t_cD1, mode=thresh_mode)
#     cD1_thresh = np.zeros(len(cD1))
#     t_cD2 = np.divide(np.median(abs(cD2)) * np.sqrt(2 * np.log(len(cD2))), 0.6745)
#     #cD2_thresh = pywt.threshold(cD2, t_cD2, mode=thresh_mode)
#     cD2_thresh = np.zeros(len(cD2))
#     #t_cD3 = np.divide(np.median(abs(cD3)) * np.sqrt(2 * np.log(len(cD3))), 0.6745)
#     #cD3_thresh = pywt.threshold(cD3, t_cD3, mode=thresh_mode)
#     #cD3_thresh = np.zeros(len(cD3))
#     cD3_thresh = cD3

#     t_cA3 = np.divide(np.median(abs(cA3)) * np.sqrt(2 * np.log(len(cA3))), 0.6745)
#     #cA3_thresh = pywt.threshold(cA3, t_cA3, mode=thresh_mode)
#     cA3_thresh = np.zeros(len(cA3))

#     # Wavelet Reconstruciton
#     coeffs_denoised = cA3, cD3_thresh, cD2_thresh, cD1_thresh
#     samples_denoised = pywt.waverec(coeffs_denoised, 'coif1', mode='zero')

#     # TK Energy Operator
#     samples_filtered_tk = []
#     for index, val in enumerate(samples_denoised[1:-1]):
#         samples_filtered_tk.append((val * val) - (samples_denoised[index - 1] * samples_denoised[index + 1]))
#     samples_filtered_tk.insert(0, samples_filtered_tk[0])
#     samples_filtered_tk.append(samples_filtered_tk[-1])

#     if plot:
#         fig_pipe = plt.figure()
#         ax10 = fig_pipe.add_subplot(511)
#         ax10.plot(ecg)
#         ax11 = fig_pipe.add_subplot(512)
#         ax11.plot(samples)
#         ax12 = fig_pipe.add_subplot(513)
#         ax12.plot(samples_stand)
#         ax13 = fig_pipe.add_subplot(514)
#         ax13.plot(samples_denoised)
#         ax14 = fig_pipe.add_subplot(515)
#         ax14.plot(samples_filtered_tk)

#         fig_wave = plt.figure()
#         ax21 = fig_wave.add_subplot(421)
#         ax21.plot(cD1)
#         ax22 = fig_wave.add_subplot(423)
#         ax22.plot(cD2)
#         ax23 = fig_wave.add_subplot(425)
#         ax23.plot(cD3)
#         ax24 = fig_wave.add_subplot(427)
#         ax24.plot(cA3)
#         ax25 = fig_wave.add_subplot(422)
#         ax25.plot(cD1_thresh)
#         ax26 = fig_wave.add_subplot(424)
#         ax26.plot(cD2_thresh)
#         ax27 = fig_wave.add_subplot(426)
#         ax27.plot(cD3_thresh)
#         ax28 = fig_wave.add_subplot(428)
#         ax28.plot(cA3_thresh)

#         plt.ion()
#         plt.show()  # block=False)
#         # plt.pause(5)
#         fig_pipe.clear()
#         fig_wave.clear()
#         plt.close('all')

#     return samples

    # # Filter for True R-Peaks
    # if len(local_peaks) == 1:  # Take it if one
    #     r_peaks.append(local_peaks[0] + i)
    # if len(local_peaks) == 2:  # Compare them if two
    #     if (local_peaks[1] - local_peaks[0]) < 75:  # 300ms in Samples
    #         if samples_window[local_peaks[1]] > samples_window[local_peaks[0]]:
    #             r_peaks.append(local_peaks[1] + i)
    #         else:
    #             r_peaks.append(local_peaks[0] + i)
    #     else:
    #         r_peaks.append(local_peaks[0] + i)
    #         r_peaks.append(local_peaks[1] + i)
    # elif len(local_peaks) > 2:  # Big Comparison
    #     mask = np.zeros(len(local_peaks))
    #     for lpc, local_peak in enumerate(local_peaks[1:-1]):
    #         if abs(local_peaks[lpc - 1] - local_peak) and abs(local_peaks[lpc + 1] - local_peak) >= 75:
    #             mask[lpc + 1] = 1  # Andere Werte sind fern von local peak
    #         else:  # Andere Werte sind nah an local peak
    #             if (samples_window[local_peaks[lpc + 2]] > samples_window[local_peak]) or (samples_window[local_peaks[lpc]] > samples_window[local_peak]):
    #                 mask[lpc + 1] = 0  # Andere Werte sind größer als local peak
    #             else:
    #                 mask[lpc + 1] = 1  # Andere Werte sind kleiner als local peak
    #     if mask[1] == 0:  # Ausnahme Erster Wert abhängig von Zweitem Wert, sonst Null lassen
    #         mask[0] = 1
    #     if mask[-2] == 0:  # Ausnahme Letzter Wert abhängig von vorletztem Wert, sonst Null lassen
    #         mask[-1] = 1
    #     local_peaks = np.array(local_peaks)
    #     mask_bool = mask.astype(bool)
    #     local_peaks = local_peaks[mask_bool]  # Update local peaks list
    #     for local_peak in local_peaks:
    #         r_peaks.append(local_peak + i)

    # up_crossings = np.argwhere(np.diff(samples_window > t_val, prepend=False))[:, 0]
    # down_crossings = np.argwhere(np.diff(samples_window > t_val, append=False))[:, 0]  # returns two vals
    # thresholds = np.diff(samples_window > t_val, prepend=False)
    # r_peaks.append(up_crossings[0] + i)
    # # True R Peak Detection
    # # take two peaks which are close to each other, relative to the the window length and calc the midposition, thats the r_peaks
    # for peak in local_peaks:
    #     # Midpoint calc
    #     r_peaks.append(ind + i)  # i equals the global pos
    #         # Condition: R-peaks has a rising and a falling edge
    # if len(local_peaks) < 2:
    #     continue

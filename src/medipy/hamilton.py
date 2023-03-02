'''
hamilton.py
'''
import scipy.signal as sc
import numpy as np

# Class with main func?


def preprocessing(samples, sample_rate):
    '''
    This method preprocesses a ecg signal for hamilton r peak detection algorithm
    '''
    # Akausaler Bandpass Filter nach [Sedghamiz2014]
    f1 = 5 / sample_rate
    f2 = 15 / sample_rate
    b, a = sc.butter(1, [f1 * 2, f2 * 2], btype='bandpass')
    samples_filtered = sc.filtfilt(b, a, samples)

    # Akusaler 5-Point-Differentiator nach [Pan1985]
    b = [1, 2, 0, -2, -1]
    b = [x * (1 / 8) * sample_rate for x in b]
    a = [1]
    samples_diff = sc.filtfilt(b, a, samples_filtered)

    # Betragsbildung (Rectifier)
    samples_rect = abs(samples_diff)

    # Akausaler ungerader 84ms Mittelwertfilter
    b = np.ones(int(0.084 * sample_rate))  # Filterbreite 21 Samples, damit ungerade
    b = b / int(0.084 * sample_rate)
    a = [1]
    samples_ma = sc.filtfilt(b, a, samples_rect)

    return samples_ma


def detect_r_peaks(samples, preprocessed, sample_rate, least_distance=0.2, th_coefficient=0.189, th_search_back=0.3):
    '''
    This method by hamilton detects all r-peaks of an ecg signal
    '''
    # Initiation
    noisy_peaks = []
    noisy_peaks_average = 0.0
    safe_peaks = []
    safe_peaks_average = 0.0
    qrs_peaks = [0]
    rr_intervals = []
    rr_intervals_average = 0.0
    threshold = 0.0
    counter = 0
    index = []
    peaks = []

    # Hamilton Algorithmus
    for counter in range(len(preprocessed)):
        # Peakdetektion
        if counter > 3 and counter < len(preprocessed) - 5:
            if preprocessed[counter - 1] < preprocessed[counter] and preprocessed[counter + 1] < preprocessed[counter]:
                peak = counter
                peaks.append(counter)
                # R4 & R1
                if preprocessed[counter] > threshold and (peak - qrs_peaks[-1]) > least_distance * sample_rate:
                    slope_neg_peak = np.diff(preprocessed[peak:peak + 2])  # 2 Samples = 8ms ANTIKAUSAL
                    slope_pos_peak = np.diff(preprocessed[peak - 2:peak])  # 2 Samples = 8ms KAUSAL
                    slope_pos_previous = np.diff(preprocessed[qrs_peaks[-1] - 2:qrs_peaks[-1]])  # 4 Samples = 8ms KAUSAL
                    if slope_pos_previous.size == 0:
                        slope_pos_previous = slope_pos_peak
                    # R2
                    if (np.mean(slope_pos_peak) > 0) and (np.mean(slope_neg_peak) < 0):
                        # Ja, dann weiter, ansonsten BSL Noise, go to R2 else
                        if (peak - qrs_peaks[-1]) < 0.360 * sample_rate and np.max(slope_pos_peak) < th_search_back * np.max(slope_pos_previous):
                            # Ja, dann T-Wave Noise, NPL Update und TH Update und Weitermachen
                            noisy_peaks.append(preprocessed[peak])
                            if len(noisy_peaks) > 8:
                                noisy_peaks.pop(0)
                            noisy_peaks_average = np.mean(noisy_peaks)  # Noisy Peak Average
                            # Noisy Peaks Average + 0.45 Parameter * (True Peaks Average - Noisy Peaks Average)
                            threshold = noisy_peaks_average + th_coefficient * (safe_peaks_average - noisy_peaks_average)
                            counter += 1
                            continue
                        # Peak ist QRS Peak
                        qrs_peaks.append(peak)
                        index.append(counter)
                        safe_peaks.append(preprocessed[peak])
                        if len(safe_peaks) > 8:
                            safe_peaks.pop(0)
                        safe_peaks_average = np.mean(safe_peaks)
                        # R5
                        if rr_intervals_average != 0.0:
                            if qrs_peaks[-1] - qrs_peaks[-2] > 1.5 * rr_intervals_average:  # Letztes RR Intervall größer als 1.5 mal RR Interval Average
                                missed_peaks = peaks[index[-2] + 1:index[-1]]  # Suche nach Peaks zwischen dem jetzigen und dem letzten
                                for missed_peak in missed_peaks:  # Gab es da einen Peak der 360ms vom letzten entfernt war? Thresshold nur die Hälfte?
                                    if missed_peak - peaks[index[-2]] > int(0.360 * sample_rate) and preprocessed[missed_peak] > 0.3 * threshold:
                                        qrs_peaks.insert(-2, missed_peak)  # Füge ihn hinzu TRUE PEAK AVERAGE? # MODIFIED INSTEAD SORT()
                                        safe_peaks.insert(-2, missed_peak)  # MODIFIED: Average neu 'CHECK THRESSHOLD
                                        if len(safe_peaks) > 8:
                                            safe_peaks.pop(0)
                                        safe_peaks_average = np.mean(safe_peaks)  # True Peak Average 'AUCH Noise peaks akutalisieren!
                                        break
                        # R5
                        if len(qrs_peaks) > 2:
                            rr_intervals.append(qrs_peaks[-1] - qrs_peaks[-2])
                            if len(rr_intervals) > 8:
                                rr_intervals.pop(0)
                            rr_intervals_average = int(np.mean(rr_intervals))  # RR Interval Average
                    # R2
                    else:
                        noisy_peaks.append(preprocessed[peak])
                        if len(noisy_peaks) > 8:
                            noisy_peaks.pop(0)
                        noisy_peaks_average = np.mean(noisy_peaks)  # Noisy Peak Average
                # R4 & R1
                else:
                    noisy_peaks.append(preprocessed[peak])
                    if len(noisy_peaks) > 8:
                        noisy_peaks.pop(0)
                    noisy_peaks_average = np.mean(noisy_peaks)  # Noisy Peak Average
                # Adjusting Threshholds
                threshold = noisy_peaks_average + th_coefficient * (safe_peaks_average - noisy_peaks_average)  # Noisy Peaks Average + 0.45 Parameter * (True Peaks Average - Noisy Peaks Average)
                counter += 1

    # Remove the inital value
    qrs_peaks.pop(0)

    return qrs_peaks


def refine_r_peaks(samples, qrs_peaks, refine_samples):
    '''
    Refines the QRS peaks by looking at the samples around the peak and taking the maximum value as the new peak.
    '''
    # Refining the peaks
    refine_area = int(refine_samples / 2)
    refined_peaks = []
    for qrs_peak in qrs_peaks:
        refined_peaks.append(np.argmax(samples[qrs_peak - refine_area:qrs_peak + refine_area]) + (qrs_peak - refine_area))
    qrs_peaks = refined_peaks

    return qrs_peaks

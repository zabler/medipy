'''
r_peak_detection.py
'''

import numpy as np


def hamilton(samples, preprocessed, sample_rate, least_distance=0.2, th_coefficient=0.189, th_search_back=0.3, refine_ms=24):
    '''
    This method by hamilton detects all r-peaks of an ecg signal
    '''
    # Initiation
    period_ms = int((1 / sample_rate) * 1000)
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

    # Initial Phase
    # Remove the inital value + 2 first detected peaks
    qrs_peaks.pop(0)
    qrs_peaks.pop(0)
    qrs_peaks.pop(0)

    # Refining
    if refine_ms is not None:
        refine_area = int((refine_ms / period_ms) / 2)
        refined_peaks = []
        for qrs_peak in qrs_peaks:
            refined_peaks.append(np.argmax(samples[qrs_peak - refine_area:qrs_peak + refine_area]) + (qrs_peak - refine_area))
        qrs_peaks = refined_peaks

    #  R Peaks an Grid anpassen
    return [qrs_peak * period_ms for qrs_peak in qrs_peaks]

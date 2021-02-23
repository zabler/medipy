'''
ecg.py
'''
import abc
import math
import pickle
import nolds
import numpy as np
from scipy import signal
from ecgdetectors import Detectors
from matplotlib import pyplot as plt
from astropy.timeseries import LombScargle


class Ecg(metaclass=abc.ABCMeta):
    '''
    This is the ecg interface class
    '''
    def __init__(self):
        # Signal Parameter
        self.samples = []
        self.preprocessed = []
        self.sample_rate = None
        self.period_ms = None
        self.grid = None
        self.r_peaks = []
        self.rr_intervals = []
        self.rr_artefacts = []
        self.unplausible_no_data = 0
        self.unplausible_not_normal = 0
        self.unplausible_artefacts = 0
        self.plausible=0

    def save_object(self, file):
        '''
        saves the ecg object as pickle
        '''
        with open(file, 'wb') as output_file:
            pickle.dump(self, output_file, pickle.HIGHEST_PROTOCOL)

    def load_object(file):
        '''
        returns an ecg object from a pickle
        '''
        with open(file, 'rb') as input_file:
            return pickle.load(input_file)
        
    def preprocessor_hamilton(self):
        '''
        preprocesses a ecg signal for hamilton algorithm
        '''
        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)
        
        # Akausaler Bandpass Filter nach [Liu2017]
        f1 = 5/self.sample_rate
        f2 = 15/ self.sample_rate
        b, a = signal.butter(1, [f1 * 2, f2 * 2], btype='bandpass')
        self.samples_filtered = signal.filtfilt(b, a, self.samples)
        
        # w, h = signal.freqz(b, a)
        # x1 = w * 250 * 1.0 / (2 * np.pi)
        # y1 = 20 * np.log10(abs(h))
        # fig, ax1 = plt.subplots(figsize=(12, 6))
        # ax1.semilogx(x1, y1, color='black')
        # ax1.axhline(y=-3, color=red)
        # ax1.axvline(x=5, color=red)
        # ax1.axvline(x=15, color=red)
        # ax1.set_ylabel('Amplitudengang [dB]')
        # ax1.set_xlabel('Frequenz [Hz]')
        # ax1.set_title('Frequenzgang Butterworth Bandpassfilter 1.Ordnung')
        # ax1.grid(which='both', linestyle='-', color='grey')
        # ax2 = ax1.twinx()
        # angles = np.unwrap(np.angle(h, deg=True))
        # ax2.semilogx(x1, angles, color=blue)
        # ax2.set_ylabel('Phasengang (Grad)', color=blue)
        # ax2.grid()
        # ax2.axis('tight')
        # plt.draw()
        # w, gd = signal.group_delay((b, a))
        # #x = w * 250 * 1.0 / (2 * np.pi)
        # fig = plt.figure(figsize=(12, 6))
        # plt.plot(w, gd)
        # plt.title('Groupdelay')
        # plt.draw()

        # Akusaler 5-Point-Differentiator nach [Pan1985]
        b = [1, 2, 0, -2, -1]
        b = [x * (1 / 8) * self.sample_rate for x in b]
        a=[1]
        self.samples_diff = signal.filtfilt(b, a, self.samples_filtered)
        
        # w, h = signal.freqz(b, a=1)
        # x1 = w * 250 * 1.0 / (2 * np.pi)
        # y1 = 20 * np.log10(abs(h))
        # plt.figure(figsize=(10,5))
        # plt.semilogx(x1, y1)
        # plt.ylabel('Amplitudengang [dB]')
        # plt.xlabel('Frequenz [Hz]')
        # plt.title('Frequenzgang akausaler 5-Punkt-Differentiator')
        # plt.grid(which='both', linestyle='-', color='grey')
        # plt.draw()

        # Betragsbildung (Rectifier)
        self.samples_rect = abs(self.samples_diff)

        # Akausaler 80ms Mittelwertfilter
        b = np.ones(int(0.084*self.sample_rate)) # Filterbreite 21 Samples, damit ungerade
        b = b/int(0.084*self.sample_rate)
        a = [1]
        samples_ma = signal.filtfilt(b, a, self.samples_rect)
        self.preprocessed = samples_ma

    def r_peak_detector_hamilton(self, least_distance=0.2, th_coefficient=0.189, th_search_back=0.3, refine_ms=80):
        '''
        This method detects all qrs beats of an ecg signal by hamilton algorithm
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
        for counter in range(len(self.preprocessed)):

            # Peakdetektion
            if counter > 3 and counter < len(self.preprocessed) - 5:
                if self.preprocessed[counter - 1] < self.preprocessed[counter] and self.preprocessed[counter + 1] < self.preprocessed[counter]:
                    peak = counter
                    peaks.append(counter)

                    # R4 & R1
                    if self.preprocessed[counter] > threshold and (peak - qrs_peaks[-1]) > least_distance * self.sample_rate:
                        slope_neg_peak = np.diff(self.preprocessed[peak:peak+2]) # 2 Samples = 8ms ANTIKAUSAL
                        slope_pos_peak = np.diff(self.preprocessed[peak - 2:peak]) # 2 Samples = 8ms KAUSAL
                        slope_pos_previous = np.diff(self.preprocessed[qrs_peaks[-1] - 2:qrs_peaks[-1]])  # 4 Samples = 8ms KAUSAL
                        if slope_pos_previous.size ==0:
                            slope_pos_previous=slope_pos_peak
                        # R2
                        if (np.mean(slope_pos_peak) > 0) and (np.mean(slope_neg_peak) < 0):
                            # Ja, dann weiter, ansonsten BSL Noise, go to R2 else
                            if (peak - qrs_peaks[-1]) < 0.360*self.sample_rate and np.max(slope_pos_peak) < 0.5*np.max(slope_pos_previous):
                                # Ja, dann T-Wave Noise, NPL Update und TH Update und Weitermachen
                                noisy_peaks.append(self.preprocessed[peak])
                                if len(noisy_peaks) > 8:
                                    noisy_peaks.pop(0)
                                noisy_peaks_average = np.mean(noisy_peaks)  # Noisy Peak Average
                                threshold = noisy_peaks_average + th_coefficient*(safe_peaks_average-noisy_peaks_average) #Noisy Peaks Average + 0.45 Parameter * (True Peaks Average - Noisy Peaks Average)
                                counter += 1
                                continue
                            # Peak ist QRS Peak
                            qrs_peaks.append(peak)
                            index.append(counter)
                            safe_peaks.append(self.preprocessed[peak])
                            if len(safe_peaks) > 8: 
                                safe_peaks.pop(0)
                            safe_peaks_average = np.mean(safe_peaks)

                            # R5
                            if rr_intervals_average != 0.0:
                                if qrs_peaks[-1] - qrs_peaks[-2] > 1.5 * rr_intervals_average:  #Letztes RR Intervall größer als 1.5 mal RR Interval Average
                                    missed_peaks = peaks[index[-2] + 1:index[-1]]  # Suche nach Peaks zwischen dem jetzigen und dem letzten
                                    for missed_peak in missed_peaks:  # Gab es da einen Peak der 360ms vom letzten entfernt war? Thresshold nur die Hälfte?
                                        if missed_peak - peaks[index[-2]] > int(0.360 * self.sample_rate) and self.preprocessed[missed_peak] > 0.3 * threshold:
                                            qrs_peaks.insert(-2, missed_peak)  # Füge ihn hinzu TRUE PEAK AVERAGE? # MODIFIED INSTEAD SORT()
                                            safe_peaks.insert(-2, missed_peak)  # MODIFIED: Average neu 'CHECK THRESSHOLD
                                            if len(safe_peaks) > 8:
                                                safe_peaks.pop(0)
                                            safe_peaks_average = np.mean(safe_peaks) # True Peak Average 'AUCH Noise peaks akutalisieren!
                                            break

                            # R5             
                            if len(qrs_peaks) > 2:
                                rr_intervals.append(qrs_peaks[-1] - qrs_peaks[-2])
                                if len(rr_intervals) > 8:
                                    rr_intervals.pop(0)
                                rr_intervals_average = int(np.mean(rr_intervals))  # RR Interval Average
                        #R2
                        else:
                            noisy_peaks.append(self.preprocessed[peak])
                            if len(noisy_peaks) > 8:
                                noisy_peaks.pop(0)
                            noisy_peaks_average = np.mean(noisy_peaks) # Noisy Peak Average
                    # R4 & R1
                    else:
                        noisy_peaks.append(self.preprocessed[peak])
                        if len(noisy_peaks) > 8:
                            noisy_peaks.pop(0)
                        noisy_peaks_average = np.mean(noisy_peaks) # Noisy Peak Average

                    #Adjusting Threshholds
                    threshold = noisy_peaks_average + th_coefficient*(safe_peaks_average-noisy_peaks_average) #Noisy Peaks Average + 0.45 Parameter * (True Peaks Average - Noisy Peaks Average)
                    counter += 1
        
        #  Initial Wert + ersten zwei Werte entfernen
        qrs_peaks.pop(0)
        qrs_peaks.pop(0)
        qrs_peaks.pop(0)

        # Refining
        if refine_ms is not None:
            refine_area = int((refine_ms/self.period_ms)/2)
            refined_peaks=[]
            for qrs_peak in qrs_peaks:
                refined_peaks.append(np.argmax(self.samples[qrs_peak - refine_area:qrs_peak+refine_area])+(qrs_peak-refine_area))
            qrs_peaks = refined_peaks

        #  R Peaks an Grid anpassen
        self.r_peaks = [qrs_peak * self.period_ms for qrs_peak in qrs_peaks]

    def rr_interval_calculator(self):
        '''
        This method calculates all rr intervals of a given r_peak list
        '''
        for k in range(1, len(self.r_peaks)):
            self.rr_intervals.append(math.ceil(self.r_peaks[k] - self.r_peaks[k - 1]))
    
    def rr_interval_kubios_artefact_detector(self, rr_intervals):
        '''
        This method checks for artefacts with an modified kubios algorithm
        '''
        alpha = 5.2
        drr_intervals = np.diff(rr_intervals)
        drr_intervals = np.insert(drr_intervals, 0, 0)
        #MODIFIED NOT 90 SUROUNDING, ALL OF 5 MIN WINDOOW FOR TH1
        threshold_1 = ((np.quantile(rr_intervals, 0.75) - np.quantile(rr_intervals, 0.25)) / 2) * alpha
        drr_intervals = drr_intervals/threshold_1
        mrr_intervals = []
        for index, rr_interval in enumerate(rr_intervals):
            if index < 5 or index > (len(rr_intervals) - 5):
                mrr_interval = rr_interval - np.median(rr_intervals)
                if mrr_interval < 0:
                    mrr_intervals.append(2 * mrr_interval)
                else:
                    mrr_intervals.append(mrr_interval)
            else:
                mrr_interval = rr_interval - np.median(rr_intervals[index - 5:index + 5])
                if mrr_interval < 0:
                    mrr_intervals.append(2 * mrr_interval)
                else:
                    mrr_intervals.append(mrr_interval)
        #MODIFIED NOT 90 SUROUNDING, ALL OF 5 MIN WINDOOW FOR TH2
        threshold_2 = ((np.quantile(mrr_intervals, 0.75) - np.quantile(mrr_intervals, 0.25)) / 2) * alpha 
        mrr_intervals = mrr_intervals / threshold_2
        
        ectopic_beats = 0
        missed_beats = 0
        extra_beats = 0
        long_short_intervals = 0
        const_1 = 0.13
        const_2 = 0.17

        for index in range(len(drr_intervals)-2):
            if abs(drr_intervals[index]) > 1:
                s11 = drr_intervals[index]
                if s11 > 0:
                    s12 = max(drr_intervals[index - 1:index + 1])
                else:
                    s12 = min(drr_intervals[index - 1:index + 1])
                if (s11 > 1 and s12 < const_1 * s11 - const_2) or (s11 < -1 and s12 > -const_1 * s11 + const_2):
                    ectopic_beats +=1
                    continue
                else:
                    if abs(drr_intervals[index]) > 1 or abs(mrr_intervals[index]) > 3:
                        if (np.sign(drr_intervals[index])*drr_intervals[index+1]<-1) or (abs(mrr_intervals[index]) > 3) or (np.sign(drr_intervals[index])*drr_intervals[index+2]<-1):
                            if abs((rr_intervals[index]/2)-mrr_intervals[index]) < threshold_2:
                                missed_beats +=1
                                continue
                            elif abs(rr_intervals[index]+rr_intervals[index]-mrr_intervals[index]) < threshold_2:
                                extra_beats +=1
                                continue
                            elif (np.sign(drr_intervals[index])*drr_intervals[index+1]<-1) or (abs(mrr_intervals[index]) > 3):
                                long_short_intervals += 1
                                continue
                            elif (np.sign(drr_intervals[index])*drr_intervals[index+2]<-1):
                                long_short_intervals += 2
                                continue
        
        # WITHOUT ECTOPICS -> NOT COUNT AS ARTEFACT IN CASE OF EPILEPSY
        self.rr_artefacts.append(ectopic_beats + missed_beats + extra_beats + long_short_intervals)
        return ectopic_beats + missed_beats + extra_beats + long_short_intervals

    def rr_plausibility_check(self, rr_intervals, window=300, güte=0.1):
        '''
        This method checks if the window has plausible values
        (1) Gibt es überhaupt Werte? Ja weiter zu (2), Nein False
        (2) Annahme: Stationarität und Normalverteilung 
        Liegt die Anzahl der detektieren R-Zacken im Bereich des theoretischen? Ja, weiter zu(3), Nein, return False
        (3) Ist der Anteil an Artefakten im 5min Fenster kleinergleich 10%?  ja, return True, Nein retrun False
        '''
        
        if len(rr_intervals) == 0:
            self.unplausible_no_data+=1
            return False
        else:
            # MEDIAN EVTL. ROBUSTER ALS MEAN, aber mean setzt verteilung vorraus, median symmetrie?
            rr_median = np.median(rr_intervals) 
            beats_theoretical = (window*1000) / rr_median
            beats_actual = len(rr_intervals) + 1
            if beats_actual not in range(int(beats_theoretical * (1 - güte)), int(beats_theoretical * (1 + güte))):
                self.unplausible_not_normal+=1
                return False
            else:
                if self.rr_interval_kubios_artefact_detector(rr_intervals) > güte * len(rr_intervals):
                    self.unplausible_artefacts+=1
                    return False
                else:
                    self.plausible +=1
                    return True

    def hrv_features_time(self, rr_intervals):
        # Statistische
        rr_median = np.median(rr_intervals)
        rr_mean = np.mean(rr_intervals)
        sdnn = np.std(rr_intervals, ddof=1)
        rmssd = np.sqrt(np.mean(np.diff(rr_intervals)** 2))
        nn50 = sum(np.abs(np.diff(rr_intervals)) > 50)
        pnn50 = (nn50 / len(rr_intervals)) * 100
        #nn20
        #pNN20
        #sdsd
        #cvsd
        #cvnn
        #madNN
        #mcvNN

        # Herzrate
        hr = np.divide(60000, rr_intervals)
        hr_mean = np.mean(hr)
        hr_max_min = max(hr)-min(hr)
        # hr_std = np.std(hr, ddof=1)
        # hr_max = max(hr)
        # hr_min = min(hr)

        #Geometrische
        tri = len(rr_intervals) / max(np.histogram(rr_intervals, bins=range(int(min(rr_intervals)), int(max(rr_intervals)), 8))[0])
        #tinn
        
        time_features = {
        'RR_MEDIAN': rr_median, 
        'RR_MEAN': rr_mean,
        'SDNN': sdnn,
        'RMSSD': rmssd,
        'NN50': nn50,
        'pNN50': nn50,
        'HR_MEAN': hr_mean,
        'HR_MAX_MIN': hr_max_min,
        'TRI': tri,
        }
        return time_features

    def hrv_features_frequency(self, rr_intervals):
        '''
        This method calculates all frequency features
        '''
        nni_tmstp = np.cumsum(rr_intervals) / 100
        timestamp_list = nni_tmstp - nni_tmstp[0]
        #why not scipy lomb
        freq, psd = LombScargle(timestamp_list, rr_intervals,
                                normalization='psd').autopower(minimum_frequency=0.04,
                                                               maximum_frequency=0.40)
        
        #vlf
        lf_indexes = np.logical_and(freq >= 0.04, freq < 0.15)
        hf_indexes = np.logical_and(freq >= 0.15, freq < 0.40)
        lf_power = np.trapz(y=psd[lf_indexes], x=freq[lf_indexes])
        hf_power = np.trapz(y = psd[hf_indexes], x = freq[hf_indexes])
        lf_peak = freq[lf_indexes].max()
        hf_peak = freq[hf_indexes].max()
        lf_hf_ratio = lf_power / hf_power
        lf_nu = (lf_power / (lf_power + hf_power)) * 100
        hf_nu = (hf_power / (lf_power + hf_power)) * 100

        freqency_features = {
        'LF_POWER': lf_power,
        'HF_POWER': hf_power,
        'LF_HF_RATIO': lf_hf_ratio,
        'LF_NU': lf_nu,
        'HF_NU': hf_nu,
        'LF_PEAK': lf_peak,
        'HF_PEAK': hf_peak
        }
    
        return freqency_features

    def hrv_features_nonlinear(self, rr_intervals):
        diff_nn_intervals = np.diff(rr_intervals)
        sd1 = np.sqrt(np.std(diff_nn_intervals, ddof=1) ** 2 * 0.5)
        sd2 = np.sqrt(2 * np.std(rr_intervals, ddof=1) ** 2 - 0.5 * np.std(diff_nn_intervals, ddof=1) ** 2)
        T = 4 * sd1
        L = 4 * sd2
        csi = sd2 / sd1
        cvi = np.log10(L * T)
        modified_csi = L ** 2 / T
        samp_en = nolds.sampen(rr_intervals, emb_dim=2)

        # df_alpha_1 nolds
        # df_alpha_2 nolds
        # cd_d2 nolds
        # ap_en
        # check neurokit

        nonlinear_features = {
            'SD1': sd1,
            'SD2': sd2,
            'CSI': csi,
            'CVI': cvi,
            'MODIFIED_CSI': modified_csi,
            'SAMPEN': samp_en
        }

        return nonlinear_features

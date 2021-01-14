'''
ecg.py
'''
import abc
import math
import pickle
import numpy as np
from scipy import signal
from ecgdetectors import Detectors
from matplotlib import pyplot as plt

class Ecg(metaclass=abc.ABCMeta):
    '''
    This is the ecg interface class
    '''
    def __init__(self):
        self.samples = []
        self.preprocessed = []
        self.sample_rate = None
        self.period_ms = None
        self.grid = None
        self.r_peaks = []
        self.rr_intervals = []
        self.rr_missings = []

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
        # Detrending
        self.samples = signal.detrend(self.samples, axis=-1, type='linear', bp=0, overwrite_data=False)

        # Bandpass Filter
        f1 = 8/self.sample_rate
        f2 = 16 / self.sample_rate
        b, a = signal.butter(1, [f1 * 2, f2 * 2], btype='bandpass')
        samples_filtered = signal.lfilter(b, a, self.samples)

        # Differentiator
        samples_diff = np.diff(samples_filtered)

        # Rectifier
        samples_rect = abs(samples_diff)

        # 80 ms Moving Average
        b = np.ones(int(0.08*self.sample_rate))
        b = b/int(0.08*self.sample_rate)
        a = [1]
        samples_ma = signal.lfilter(b, a, samples_rect)
        samples_ma[0:len(b) * 2] = 0  # Einschwingtiefe
        self.preprocessed = samples_ma

    def r_peak_detector_hamilton(self, least_distance=0.2, th_coefficient=0.45):
        '''
        This method detects all qrs patterns of ecg signal by hamilton algorithm
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
        for counter in range(len(self.preprocessed)): # evtl problematisch da kein 250 Hz grid fenster ENUMERATE

            # Peakdetektion
            if counter > 0 and counter < len(self.preprocessed) - 1:
                if self.preprocessed[counter - 1] < self.preprocessed[counter] and self.preprocessed[counter + 1] < self.preprocessed[counter]:
                    peak = counter
                    peaks.append(counter)

                    # R4 & R1
                    if self.preprocessed[counter] > threshold and (peak - qrs_peaks[-1]) > least_distance * self.sample_rate: # MODIFIED to least_distance = 0.195
                        qrs_peaks.append(peak)
                        index.append(counter)
                        safe_peaks.append(self.preprocessed[peak])
                        if len(safe_peaks) > 8:  #MODIFIED n to s
                            safe_peaks.pop(0)
                        safe_peaks_average = np.mean(safe_peaks)

                        # R5
                        if rr_intervals_average != 0.0:
                            if qrs_peaks[-1] - qrs_peaks[-2] > 1.5 * rr_intervals_average:  #Letztes RR Intervall größer als 1.5 mal RR Interval Average
                                missed_peaks = peaks[index[-2] + 1:index[-1]]  # Suche nach Peaks zwischen dem jetzigen und dem letzten
                                for missed_peak in missed_peaks:  # Gab es da einen Peak der 360ms vom letzten entfernt war? Thresshold nur die Hälft?
                                    if missed_peak - peaks[index[-2]] > int(0.360 * self.sample_rate) and self.preprocessed[missed_peak] > 0.5 * threshold:
                                        qrs_peaks.insert(-2, missed_peak)  # Füge ihn hinzu TRUE PEAK AVERAGE? # MODIFIED INSTEAD SORT()
                                        safe_peaks.insert(-2, missed_peak)  # MODIFIED: Average neu 'CHECK THRESSHOLD
                                        if len(safe_peaks) > 8:
                                            safe_peaks.pop(0)
                                        safe_peaks_average = np.mean(safe_peaks) # True Peak Average
                                        break

                        # R5             
                        if len(qrs_peaks) > 2:
                            rr_intervals.append(qrs_peaks[-1] - qrs_peaks[-2])
                            if len(rr_intervals) > 8:
                                rr_intervals.pop(0)
                            rr_intervals_average = int(np.mean(rr_intervals))  # RR Interval Average

                    # R4
                    else:
                        noisy_peaks.append(self.preprocessed[peak])
                        if len(noisy_peaks) > 8:
                            noisy_peaks.pop(0)
                        noisy_peaks_average = np.mean(noisy_peaks) # Noisy Peak Average

                    threshold = noisy_peaks_average + th_coefficient*(safe_peaks_average-noisy_peaks_average) #Noisy Peaks Average + 0.45 Parameter * (True Peaks Average - Noisy Peaks Average)
                    counter += 1
        
        #  Initial Wert + ersten zwei Werte entfernen
        qrs_peaks.pop(0)
        qrs_peaks.pop(0)
        qrs_peaks.pop(0)

        #  R Peaks an Grid anpassen
        self.r_peaks = [qrs_peak * self.period_ms for qrs_peak in qrs_peaks]

    def rr_interval_calculator(self):
        '''
        This method calculates all rr intervals of a given r_peak list
        '''
        for k in range(1, len(self.r_peaks)):
            self.rr_intervals.append(math.ceil(self.r_peaks[k] - self.r_peaks[k - 1]))
    
    def rr_interval_wrong_data_detector(self):
        '''
        This method detects all wrong rr_intervals
        '''

    def rr_interval_wrong_data_remover(self):
        '''
        This method removes all the wrong r peaks of the wrong rr_intervals
        '''

    def rr_interval_missing_data_detector(self):
        '''
        This method detects all missing rr_intervals
        '''

        for index, rr_interval in enumerate(self.rr_intervals):
            if index <= 8:
                self.rr_missings.append(0)
                continue
            rr_threshold = np.mean(self.rr_intervals[index - 9:index]) + np.std(self.rr_intervals[index - 9:index])
            if rr_interval >= rr_threshold:
                rr_missing = int(round(rr_interval / rr_threshold)) - 1
                self.rr_missings.append(rr_missing)
            else:
                self.rr_missings.append(0)
    
            # rr bigger als umgebung
            # rr größer als mean der letzten und kommenden 8 werte
            # dann sag artefact:
            # wert davor, wert danach, daraus anzahl der benötigten R zacken berechnen
        #return liste genauso lang wie rr_interval, aber mit int angabe der Werte die vom letzen R wert bis jetzt zu interpolieren sind 
    
    def rr_interval_interpolator(self, rr_list, method=None):  #time, rr_missing,
        '''
        This method interpolates rr_intevals based on time information
        '''
        if method is None:
            return rr_list
        elif method == 'linear':
                #rr_interpolated = signal.interpolate.interp1d(,,kind=method) #kompliziert
            rr_interpolated = rr_list[~np.isnan(rr_list)]
            return rr_interpolated

    def hrv_feature_sdnn(self, rr_intervals):
        '''
        This method calculates sdnn based on a ideal rr_interval list
        '''

    def hrv_feature_sdann(self, rr_intervals):
        '''
        This method calculates sdann based on a ideal rr_interval list
        '''
        feature = np.mean(rr_intervals) 
        return feature
        #returnd value
        # to be interpoalted = missing_data_check(rr[])
        # kann evtl auch außerhalb der Funktion im calculator gefragt werdeb
        # elif summe der rr_missing aus diesem bereich == 10% aller rr_werte aus diesem Bereich
        # return np.nan als wert
        # elif
        # interpolate nach methode 
        # reutrn dann feature 
        # was wenn erster RR Wert rückwirkend interpolieren muss?, dann nein

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
        self.sample_rate = None
        self.delta_ms = None
        self.r_peaks = []
        self.rr_intervals = []
        self.rr_missing = []
    
    @abc.abstractmethod
    def data_reader(self, file):
        '''
        gets the data from files
        '''
        raise NotImplementedError
    
    @abc.abstractmethod
    def save_to_table(self, file):
        '''
        saves the data in a table format
        '''
        raise NotImplementedError

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
    
    
    def preprocessing_hamilton(self):
        '''
        preprocessing for hamilton algorithm
        '''
        # Detrending
        self.samples = signal.detrend(self.samples, axis= - 1, type='linear', bp=0, overwrite_data=False)

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
        self.preprocessed_hamilton = samples_ma
    
    def r_peak_detection_hamilton(self, least_distance=0.2, th_coefficient=0.45):
        
        # Initiation
        noisy_peaks = []
        noisy_peaks_average = 0.0
        safe_peaks = []
        safe_peaks_average = 0.0
        qrs_peaks = [0]
        rr_intervals = []
        rr_intervals_average = 0.0
        threshold = 0.0
        counter=0
        index = []
        peaks = []  

        # Hamilton Algorithmus
        for counter in range(len(self.preprocessed_hamilton)):

            # Peakdetektion
            if counter>0 and counter<len(self.preprocessed_hamilton)-1: 
                if self.preprocessed_hamilton[counter-1]<self.preprocessed_hamilton[counter] and self.preprocessed_hamilton[counter+1]<self.preprocessed_hamilton[counter]:
                    peak = counter
                    peaks.append(counter)

                    # R4 & R1
                    if self.preprocessed_hamilton[counter] > threshold and (peak-qrs_peaks[-1])>least_distance*self.sample_rate: # MODIFIED to least_distance = 0.195
                        qrs_peaks.append(peak)
                        index.append(counter)
                        safe_peaks.append(self.preprocessed_hamilton[peak])
                        if len(safe_peaks)>8: #MODIFIED n to s
                            safe_peaks.pop(0)
                        safe_peaks_average = np.mean(safe_peaks)

                        # R5 
                        if rr_intervals_average != 0.0: 
                            if qrs_peaks[-1]-qrs_peaks[-2] > 1.5*rr_intervals_average: #Letztes RR Intervall größer als 1.5 mal RR Interval Average
                                missed_peaks = peaks[index[-2]+1:index[-1]] # Suche nach Peaks zwischen dem jetzigen und dem letzten 
                                for missed_peak in missed_peaks: # Gab es da einen Peak der 360ms vom letzten entfernt war? Thresshold nur die Hälft?
                                    if missed_peak-peaks[index[-2]]>int(0.360*self.sample_rate) and self.preprocessed_hamilton[missed_peak]>0.5*threshold:
                                        qrs_peaks.insert(-2, missed_peak)  # Füge ihn hinzu TRUE PEAK AVERAGE? # MODIFIED INSTEAD SORT()
                                        safe_peaks.insert(-2, missed_peak) # MODIFIED: Average neu 'CHECK THRESSHOLD
                                        if len(safe_peaks)>8:
                                            safe_peaks.pop(0)
                                        safe_peaks_average = np.mean(safe_peaks) # True Peak Average
                                        break

                        # R5                 
                        if len(qrs_peaks)>2:
                            rr_intervals.append(qrs_peaks[-1]-qrs_peaks[-2])
                            if len(rr_intervals)>8:
                                rr_intervals.pop(0)
                            rr_intervals_average = int(np.mean(rr_intervals)) # RR Interval Average

                    # R4 
                    else: 
                        noisy_peaks.append(self.preprocessed_hamilton[peak])
                        if len(noisy_peaks)>8:
                            noisy_peaks.pop(0)
                        noisy_peaks_average = np.mean(noisy_peaks) # Noisy Peak Average

                    threshold = noisy_peaks_average + th_coefficient*(safe_peaks_average-noisy_peaks_average) #Threshold Berechnung: Noisy Peaks Average + 0.45 Parameter * (True Peaks Average - Noisy Peaks Average)

                    counter+=1
        
        #  Initial Wert + ersten zwei Werte entfernen
        qrs_peaks.pop(0)  
        qrs_peaks.pop(0)
        qrs_peaks.pop(0)
        self.r_peaks = qrs_peaks
        # plt.figure()
        # plt.plot(ecg.samples)
        # plt.plot(QRS, ecg.samples[QRS],'x')

    
    def rr_interval_calculator(self):
        '''
        This method calculates all rr intervals of a given r_peak list
        '''
        for k in range(1, len(self.r_peaks)):
            self.rr_intervals.append(math.ceil((self.r_peaks[k] - self.r_peaks[k - 1]) * self.delta_ms))
    
    def rr_interval_missing_data_detector(self):
        '''
        This method detects all rr abornamlties of a given rr list 
        and returns a list about the values which have to be interpolated 
        '''
        for index, rr in enumerate(self.rr):
            # rr bigger als umgebung
            # rr größer als mean der letzten und kommenden 8 werte
            # dann sag artefact:
            # wert davor, wert danach, daraus anzahl der benötigten R zacken berechnen
        #return liste genauso lang wie rr_interval, aber mit int angabe der Werte die vom letzen R wert bis jetzt zu interpolieren sind 
        
        pass
    
    def missing_data_check(self, rr_list, rr_missing, güte=0.1):
        '''
        This method checks if enough actual values are here
        '''
        if sum(rr_missing)/len(rr_list) <= güte:
            return True
        else:
            return False

        #return summe für rr bereich über wie vieele R-zacken interpoliert werden müssen

    def rr_interpolation(self, time, rr_list, rr_missing, method=None)
        '''
        Interpoliert rr_interval bei geggebenen rr_missung nach basierender methode
        '''
        if method == None:
            return rr_list
        elif method = linear:
                #scipy linare interpilation
            return rr_interpolated
    
    def hrv_feature_sdann(self, rr_good)
        '''
        Berechnet das entsprechende feature, nimmt eine perfekte RR_liste an, bekomme keine weiteren Informatuionen
        '''
        #returnd value
        # to be interpoalted = missing_data_check(rr[])
        # kann evtl auch außerhalb der Funktion im calculator gefragt werdeb
        # elif summe der rr_missing aus diesem bereich == 10% aller rr_werte aus diesem Bereich
        # return np.nan als wert
        # elif
        #   interpolate nach methode 
        #   reutrn dann feature 
        # was wenn erster RR Wert rückwirkend interpolieren muss?, dann nein

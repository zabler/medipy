'''
This is the ecg.py file
'''
import math
import numpy as np
from wfdb import processing
from pyedflib import highlevel
from medipy.interfaces.signal import Signal

# Physikalische Umrechnung
# Allgemeine Add Channel Funktion

class Ecg(Signal):
    '''
    This is the Ecg class
    '''

    def __init__(self):
        self.samples = []
        self.time = []
        self.tags = []
        self.sample_rate = None
        self.delta_ms = None
        self.r_peaks = []
        self.rr_intervals = []

    def from_edfplus_uklfr(self, path):
        '''
        This method reads a ecg signal from edfplus in the format of uniklinik freiburg
        '''
        # Specific Input (time effiency)
        read_list = ['EKG 1', 'EKG 2']

        # Read In
        signals, signal_headers, header = highlevel.read_edf(path, ch_names=read_list)
        annotations = header['annotations']  #format Seconds.Milliseconds / Duration in sec / Text= R-Wave or Seizure_TAG

        # Signal
        self.samples = signals[1] - signals[0]

        # Sample Rate
        self.sample_rate = signal_headers[0]['sample_rate']
        self.delta_ms = int((1/self.sample_rate)*1000)

        # Time
        self.time = list(np.arange(0,len(self.samples) * self.delta_ms, self.delta_ms,dtype='int32'))

        # Rpeaks and Tags #Hier bereits umrechnen! Seconds *1000 + Miliseconds
        def annotation_to_ms(timestring):
            sec, msec = str(timestring).split('.')
            if len(msec) == 1:
                msec = int(msec) * 100
            elif len(msec) == 2:
                msec = int(msec) * 10
            else:
                msec = int(msec)                
            return (int(sec)*1000)+msec

        tags = ['SEIZURE', 'EEGSTART', 'EEGSTOP']
        for annotation in annotations:
            if 'R-wave' in annotation[2]:
                self.r_peaks.append(annotation_to_ms(annotation[0]))
            elif any(tag in annotation[2] for tag in tags):
                self.tags.append([annotation_to_ms(annotation[0]), annotation[2]]) #implement here the position
        
        # RR Intervals
        for r_peak in range(1, len(self.r_peaks)):
            self.rr_intervals.append(math.ceil(self.r_peaks[r_peak] - self.r_peaks[r_peak - 1]))

    def resample_to_fs(self, new_fs):
        '''
        This method resamples the origin signal to a new target sampling frequency, e.g. 1000Hz
        '''
        self.samples, locations = processing.resample_sig(self.samples, self.sample_rate, new_fs)
        self.sample_rate = new_fs

    def get_special_tags(self, special_tag_string):
        '''
        metode in Signal einbetten
        '''
        special_tags = []
        for tag in self.tags:
            if special_tag_string in tag[1]:
                special_tags.append(tag)
        return special_tags

    def plot_signal(self):
        '''
        This method plots a small part of an ecg signal
        '''

    def r_peak_detector(self, detector):
        '''
        This method detects all r peaks in an ecg signal with the prefered detector
        '''

    def hrv_analysis(self):
        '''
        This method analysis the r peaks of an ecg signal
        '''

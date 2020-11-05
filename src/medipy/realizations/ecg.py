'''
This is the ecg.py file
'''
import math
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
        self.tags = []
        self.sample_rate = None
        self.r_peaks = []
        self.rr_intervals = []

    def from_edfplus_uklfr(self, path):
        '''
        This method reads a ecg signal from edfplus in the format of uni freiburg
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

        # Rpeaks and Tags
        for annotation in annotations:
            if 'R-wave' in annotation[2]:
                self.r_peaks.append(annotation[0])
            elif 'Seizure' in annotation[2]:
                self.tags.append(annotation[0])
        
        # RR Intervals
        for r_peak in range(1, len(self.r_peaks)):
            self.rr_intervals.append(math.ceil((self.r_peaks[r_peak] - self.r_peaks[r_peak - 1])*1000))

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

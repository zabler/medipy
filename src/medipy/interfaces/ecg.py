'''
ecg.py
'''
import abc
import pickle
import rdetection as rd

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

    def save_object(self, file):
        '''
        saves the ecg object as pickle
        '''
        with open(file, 'wb') as output_file:  # Overwrites any existing file.
            pickle.dump(self, output_file, pickle.HIGHEST_PROTOCOL)

    def load_object(file):
        '''
        returns an ecg object from a pickle
        '''
        with open(file, 'rb') as input_file:
            return pickle.load(input_file)
    
    def r_peak_detector(self):
        '''
        This method detects all r peaks in an ecg signal with the prefered detector
        rdetectors
        pyecg detectors
        wfdb
        paper
        '''
        self.r_peaks, self.rr_intervals = rd.Skipi(self.samples, self.sample_rate).detect()
    
    def r_peak_optimizer(self):
        # Verteilung checken
        # Outlier Removing
        # Median Plot
        pass
    
    def hrv_analyzer(self):
        pass

    @abc.abstractmethod
    def plot_signal(self):
        '''
        plots a part of the signal in an extra plot
        '''
        raise NotImplementedError
    
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
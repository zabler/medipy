'''
ecg.py
'''
import math
import pickle
import csv
#use dill or feather instead
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pyedflib import highlevel
from wfdb import processing
from ecgdetectors import Detectors
from medipy.interfaces.signal import Signal
import rdetection as rd

class Ecg(Signal):
    '''
    This is the Ecg class
    '''

    def __init__(self):
        self.samples = []
        self.tags = []
        self.sample_rate = None
        self.delta_ms = None
        self.r_peaks = []
        self.rr = []
        self.meta = None
    
    def plot_signal(self):
        '''
        plots a part of the signal in an extra plot
        
        '''
        fig = plt.figure(figsize=(10, 3))

        # Plot Data
        plt.plot(self.samples, label='EKG', linewidth=0.7, color='black')

        #Plot R-Peaks
        plt.plot(self.r_peaks, self.samples[self.r_peaks], 'x', color='red', label='Detektierte R-Zacken')

        # Plot Settings
        plt.xlabel(f'Zeit[s] @ {self.sample_rate}Hz ', fontname="Arial")
        plt.xlim(10*self.sample_rate, 20 * self.sample_rate)
        plt.ylim(-3,3) 
        plt.ylabel('EKG Einthoven II [mV]', fontname="Arial")
        plt.grid(b=True, which='major', axis='both')
        plt.legend(fontsize='x-small', bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4)
        plt.gca().set_xticks(np.arange(10 * self.sample_rate, 20 * self.sample_rate + 1, 2 * self.sample_rate))
        newtime = ['10', '12', '14', '16', '18', '20']
        plt.gca().set_xticklabels(newtime)
        plt.show(block=false)


    def from_edfplus_uklfr(self, path):
        '''
        This method reads a ecg signal from edfplus in the format of uniklinik freiburg
        '''
        # Specific Input
        read_list = ['EKG 1', 'EKG 2']

        # Read In
        signals, signal_headers, header = highlevel.read_edf(path, ch_names=read_list)
        annotations = header['annotations']  #format Seconds.Milliseconds / Duration in sec / Tag: SEIZURE_P##_S##

        # Signal
        self.samples = signals[1] - signals[0]

        # Sample Rate
        self.sample_rate = signal_headers[0]['sample_rate']
        self.delta_ms = int((1/self.sample_rate)*1000)

        # Tags
        def annotation_to_grid(self, timestring):
            sec, msec = str(timestring).split('.')
            if len(msec) == 1:
                msec = int(msec) * 100
            elif len(msec) == 2:
                msec = int(msec) * 10
            else:
                msec = int(msec)                
            return int(((int(sec) * 1000) + msec)/self.delta_ms)

        tags = ['SEIZURE']
        for annotation in annotations:
            if any(tag in annotation[2] for tag in tags):
                self.tags.append([annotation_to_grid(self, annotation[0]), annotation[2]]) 

    def r_peak_detector(self):
        '''
        This method detects all r peaks in an ecg signal with the prefered detector
        rdetectors
        pyecg detectors
        wfdb
        paper
        '''
        self.r_peaks, self.rr = rd.Skipi(self.samples, self.sample_rate).detect()
    
    def r_peak_optimizer(self):
        # Verteilung checken
        # Outlier Removing
        # Median Plot
        
        pass
        
    def save_signal(self, file):
        '''
        saves the signal object as pickle
        '''
        with open(file, 'wb') as output:  # Overwrites any existing file.
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
    
    def load_signal(file):
        '''
        loads the signal object from a pickle
        '''
        with open(file, 'rb') as input:
            return pickle.load(input)
    

    def save_to_table(self, file, meta_path):
        '''
        saves the data as a pandas table in pickle
        '''

        # Tags in Dataframe
        tags_values = []
        tags_index = []
        for tag in self.tags:
            tags_values.append(tag[1])
            tags_index.append(tag[0])
        tags_df = pd.DataFrame({'TAGS': tags_values}, index=tags_index)

        # Metas in Dataframe
        meta_df = pd.read_csv(meta_path, delimiter=';')
        meta_df = meta_df[meta_df['TAG'].isin(tags_values)].reset_index(drop=True)
        meta_df = meta_df.set_index(tags_df.index)
        self.meta = meta_df.values.tolist()

        # Alle Spalten Droppen die für weitere Verarbeitung irrelavant sind
        meta_df = meta_df.drop(columns=['MYOCLONIC','EPILEPTIC','TYP I','TYP II','TYP III','SEITE','TAG','STUDY_ID','ORIGIN_START','EMG','DELRE','DELLI','QUADRE','QUADLI','VIDEO','VIGILANZ','PERCEPTION','IMPULSE'])  
        
        # Samples in Dataframe
        samples_df = pd.DataFrame({'SAMPLES': self.samples})

        # R-Peaks in Dataframe
        rr = self.rr.insert(0, np.nan)
        r_peaks_df = pd.DataFrame({'R_PEAKS': rr}, index=self.r_peaks)

        # Merge Dataframe
        df = pd.concat([samples_df, tags_df, meta_df, r_peaks_df], ignore_index=False, axis=1)

        #Save Dataframe to pickle
        df.to_pickle(file)


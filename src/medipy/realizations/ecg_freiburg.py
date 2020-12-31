'''
ecg_freiburg.py
'''

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pyedflib import highlevel
from medipy.interfaces.ecg import Ecg

class EcgFreiburg(Ecg):
    '''
    This is the ecg freiburg realization class
    '''

    def __init__(self):
        super().__init__()
        self.tags = []
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
        plt.show()#block=False)

    def data_reader(self, path):
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
        rr_intervals = self.rr_intervals.insert(0, np.nan)
        r_peaks_df = pd.DataFrame({'R_PEAKS': rr_intervals}, index=self.r_peaks)

        # Merge Dataframe
        all_df = pd.concat([samples_df, tags_df, meta_df, r_peaks_df], ignore_index=False, axis=1)

        #Save Dataframe to pickle
        all_df.to_pickle(file)


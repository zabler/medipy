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
        self.seizures = []
        self.meta = None
        self.feature_df = None
    
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

        # Period_ms
        self.period_ms = int((1 / self.sample_rate) * 1000)
        
        # Grid
        self.grid = np.arange(0, len(self.samples) * self.period_ms, self.period_ms, dtype=int)

        # Tags
        def annotation_to_grid(self, timestring):
            sec, msec = str(timestring).split('.')
            if len(msec) == 1:
                msec = int(msec) * 100
            elif len(msec) == 2:
                msec = int(msec) * 10
            else:
                msec = int(msec)                
            return int(((int(sec) * 1000) + msec) / self.period_ms)
        
        tags = ['SEIZURE']
        for annotation in annotations:
            if any(tag in annotation[2] for tag in tags):
                self.seizures.append([annotation_to_grid(self, annotation[0]), annotation[2]])
            else:
                self.tags.append([annotation_to_grid(self, annotation[0]), annotation[2]])

        # tags = ['SEIZURE']
        # for annotation in annotations:
        #     self.tags.append([annotation_to_grid(self, annotation[0]), annotation[2]])
        #     if any(tag in annotation[2] for tag in tags):
        #         self.seizures.append([annotation_to_grid(self, annotation[0]), annotation[2]])

        # tags = ['SEIZURE']
        # for annotation in annotations:
        #     if any(tag in annotation[2] for tag in tags):
        #         self.tags.append([annotation_to_grid(self, annotation[0]), annotation[2]])
    
    def short_time_hrv_extractor(self, window=300, overlap=0, include_meta=None, export=None,):
        '''
        st_hrv_calculator renntü ber eine PAndas table
        Berechnet HRV Feature XY, basierend auf Zeitbasis als Index mit zentralem Speicherzeitpunkt, RRList im Zeitraum, und markierter Fehler in einem Zeitfenster von 5min überlappend ?
        je nach FEature und je nach markierten Fehler 
        in ecg freiburg
        '''
        # Feature Names
        time_features = ['SDNN', 'SDRR', 'SDANN', 'SDNN Index', 'RMSSD', 'NN50', 'pNN50', 'DeltaHRMaxHRMin', 'TRI', 'TINN']
        freq_features = ['VLF_Power', 'LF_Power', 'LF_rel', 'LF_Nu', 'LF_Peak', 'HF_Power', 'HF_rel', 'HF_Nu', 'HF_Peak', 'LF/HF']
        non_linear_features = ['S', 'SD1', 'SD2', 'SD1/SD2', 'ApEn', 'SampEn', 'DF_alpha1', 'DF_alpha2', 'D_2']
        
        # Samples in Dataframe
        samples_df = pd.DataFrame({'SAMPLES': self.samples}, index=self.grid)

        # R-Peaks, RR_Intervals and RR_Missings to Dataframe
        rr_intervals = self.rr_intervals
        rr_intervals = np.insert(rr_intervals, 0, int(np.mean(self.rr_intervals[0:10])))
        rr_missings = self.rr_missings
        rr_missings = np.insert(rr_missings, 0, 0)
        r_peaks_df = pd.DataFrame({'R_PEAKS': self.r_peaks,'RR_INTERVALS': rr_intervals, 'RR_MISSINGS': rr_missings}, index=self.r_peaks)
        
        # Combine in Signal Dataframe
        signal_df = pd.concat([samples_df, r_peaks_df], ignore_index=False, axis=1).fillna(value=np.nan)
        
        # Get from Signal Dataframe including time information
        rr_intervals_frame = signal_df['RR_INTERVALS'].to_numpy()
        rr_missings_frame = signal_df['RR_MISSINGS'].to_numpy()
        
        # Create Feature Dataframe
        self.feature_df = pd.DataFrame(np.nan, index=self.grid, columns=time_features+freq_features+non_linear_features)

        # Calculation WRONG
        half_window_ms = int(window / 2) * 1000
        step = int(((1-overlap)*half_window_ms*2)/self.period_ms)
        if step == 0:
            step = 1
        for time in self.grid[::step]: #Overlap damit nicht für jedes Sample ein Feature berechnet wird
            if time < half_window_ms:
                continue
            rr_intervals_window = rr_intervals_frame[time - half_window_ms:time + half_window_ms]
            rr_missings_window = rr_missings_frame[time - half_window_ms:time + half_window_ms]
            if self.missing_data_check(rr_intervals_window[~np.isnan(rr_intervals_window)], rr_missings_window[~np.isnan(rr_missings_window)], güte=0.1):
                self.feature_df.at[time, 'SDANN'] = self.hrv_feature_sdann(self.rr_interval_interpolator(rr_intervals_window, method='linear'))
        
        if include_meta is not None:

            # Tags in Dataframe
            tags_values = []
            tags_index = []
            for tag in self.tags:
                tags_values.append(tag[1])
                tags_index.append(tag[0])
            tags_df = pd.DataFrame({'TAGS': tags_values}, index=tags_index)
            # INDEX FÜR SYSTEMATISCHE OTHER SEIZURES FEHLT, MÜSSTE ANNOTIERT WERDEN
            # SIND ABER IN TAGS ENTHALTEN 
            
            # Seizures in Dataframe
            seizures_values = []
            seizures_index = []
            for seizure in self.seizures:
                seizures_values.append(seizure[1])
                seizures_index.append(seizure[0])
            seizures_df = pd.DataFrame({'SEIZURES': seizures_values}, index=seizures_index)

            # Metas of Seizures in Dataframe
            meta_df = pd.read_csv(include_meta, delimiter=';')
            meta_df = meta_df[meta_df['TAG'].isin(seizures_values)].reset_index(drop=True)
            meta_df = meta_df.set_index(seizures_df.index)
            self.meta = meta_df.values.tolist()

             # Alle Spalten Droppen die für weitere Verarbeitung irrelavant sind
            meta_df = meta_df.drop(columns=['MYOCLONIC','EPILEPTIC','TYP I','TYP II','TYP III','SEITE','TAG','STUDY_ID','ORIGIN_START','EMG','DELRE','DELLI','QUADRE','QUADLI','VIDEO','VIGILANZ','PERCEPTION','IMPULSE'])

            # Merge Dataframe
            self.feature_df = pd.concat([self.feature_df, seizures_df, meta_df], ignore_index=False, axis=1) #tags_df
        
        if export is not None:
            #Save Dataframe to pickle (evtl. zu Testzwecken als csv)
            self.feature_df.to_pickle(export)

    def missing_data_check(self, rr_list, rr_missing, güte=0.1):
        '''
        This method checks if there are enough non interpolated values
        '''
        #check 5min mindestwert an rr listen, was wenn rr_missing werte erst im anschluss kommen , plausbiel check mean rr -> Anzahl der mind peaks
        if len(rr_list) > 150: #30 BPM 5 min  
            if sum(rr_missing)/len(rr_list) <= güte:
                return True
        else:
            return False

    def plot_tag_in_ecg(self, tag, sec_pre=20, sec_post=40):
        '''
        plots a part of the signal in an extra plot
        
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Data
        plt.plot(self.samples, label='EKG', linewidth=1.0, color='black')

        #Plot R-Peaks
        plt.plot(self.r_peaks, self.samples[self.r_peaks], 'x', color='red', label='Detektierte R-Zacken')
        
        #Plot Tag
        plt.plot(0, 0, 'b--', label=f'{tag[1]}')
        plt.axvline(x=tag[0], color='b', linestyle='--')

        # Plot Bereich
        lower_limit = tag[0] - sec_pre*self.sample_rate
        upper_limit = tag[0] + sec_post * self.sample_rate
        
        # Plot Settings
        plt.xlabel(f'Zeit[s] @ {self.sample_rate}Hz ', fontname="Arial")
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-3, 3) 
        plt.ylabel('EKG Einthoven II [mV]', fontname="Arial")
        plt.grid(b=True, which='major', axis='both')
        plt.legend(fontsize='x-small', bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4)
        new_xticks = np.arange(lower_limit, upper_limit + 1, 4 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(-1 * sec_pre, sec_post, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.show()#block=False)

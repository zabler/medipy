'''
ecg_freiburg.py
'''

import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
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
        self.patient_id = None
    
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
            time = round(((int(sec) * 1000) + msec) / self.period_ms)
            if time in self.grid:
                return time
            else:
                index = np.abs(self.grid-time).argmin() 
                return self.grid[index]
        
        tags = ['SEIZURE']
        for annotation in annotations:
            if any(tag in annotation[2] for tag in tags):
                self.seizures.append([annotation_to_grid(self, annotation[0]), annotation[2]])
        
        #Patient ID
        filename = path.split('/')[-1]
        self.patient_id = filename.split('_')[0]
            
    def data_cutter(self, duration_s=None):
        STOP = duration_s * self.sample_rate
        self.samples = np.array(self.samples[0:STOP])
        self.grid = np.arange(0, duration_s * self.sample_rate * self.period_ms, self.period_ms, dtype=int)
        self.seizures = [seizure for seizure in self.seizures if seizure[0] in range(self.grid[0], self.grid[-1])]
    
    def short_term_hrv_extractor(self, window=300, overlap=0.5, include_meta=None, export=None):
        '''
        Short Term HRV Extraktor
        Iteration mittels geleitendem (Overlap=1) Fenster / oder springendem (Overlap<1) Fenster über Grid (entspricht Index des Dataframes)
        Für jedes Fenster werden die darin befindenden RR Intervalle auf Plausibilität und Fehlende Werte geprüft
        - Wenn nicht plausibel, wird nichts berechnet: np.nan
        - Wenn plausibel,
                - werden für den zentralen Gridwert (Window/2) verschiedene HRV Feature berechnet
                - zuvor wird je nach Feature und fehlender Datenmenge auf verschiedene Weise interpoliert
                - im Dataframe gespeichert
        - Anschließend wird der Dataframe mit den Signal, Tags und Meta Daten zusammengeführt, wenn Meta Daten zur Verfügung stehen
        - Der DataFrame wird unter "export" exportiert
        '''
        # Feature Names
        time_feature_names = ['RR_MEDIAN', 'RR_MEAN', 'SDNN', 'RMSSD', 'NN50', 'pNN50', 'HR_MAX_MIN', 'HR_MEAN', 'TRI']
        frequency_feature_names = ['LF_POWER', 'LF_NU', 'LF_PEAK', 'HF_POWER', 'HF_NU', 'HF_PEAK', 'LF_HF_RATIO']
        nonlinear_feature_names = ['SD1', 'SD2', 'CSI', 'CVI', 'MODIFIED_CSI', 'SAMPEN']

        # Samples in Dataframe
        samples_df = pd.DataFrame({'SAMPLES': self.samples}, index=self.grid)

        # R-Peaks, RR_Intervals and RR_Missings to Dataframe
        rr_intervals = self.rr_intervals
        rr_intervals = np.insert(rr_intervals, 0, int(np.mean(self.rr_intervals[0:10])))
        r_peaks_df = pd.DataFrame({'R_PEAKS': np.array(self.r_peaks), 'RR_INTERVALS': rr_intervals}, index=np.array(self.r_peaks))
        
        # Combine in Signal Dataframe
        signal_df = pd.concat([samples_df, r_peaks_df], ignore_index=False, axis=1).fillna(value=np.nan)
        
        # Get from Signal Dataframe including time information
        rr_intervals_frame = signal_df['RR_INTERVALS'].to_numpy()
        
        # Create Feature Dataframe
        self.feature_df = pd.DataFrame(np.nan, index=self.grid, columns=time_feature_names+frequency_feature_names+nonlinear_feature_names)

        # HRV Calculation
        half_window = int((window / 2) * self.sample_rate)
        step_size = int(math.ceil((1 - overlap) * half_window * 2))
        if step_size == 0:
            step_size = 1
        steps = np.arange(0, len(self.grid) + 1, step=step_size)
        for step in steps:
            if step < half_window or step+half_window > steps[-1]:
                continue
            rr_intervals_window = rr_intervals_frame[step - half_window:step + half_window]
            rr_intervals_window = rr_intervals_window[~np.isnan(rr_intervals_window)]
            if self.rr_plausibility_check(rr_intervals_window, window=300, güte=0.1):
                time_features = self.hrv_features_time(rr_intervals_window)
                frequency_features = self.hrv_features_frequency(rr_intervals_window)
                nonlinear_features = self.hrv_features_nonlinear(rr_intervals_window)
                features = {**time_features, **frequency_features, **nonlinear_features}
                for feature in features:
                    self.feature_df.at[self.grid[step], feature] = features[feature]
                
        # Append Meta Infos
        if include_meta is not None:

            # Seizures in Dataframe
            seizures_values = []
            seizures_index = []
            for seizure in self.seizures:
                seizures_values.append(seizure[1])
                seizures_index.append(seizure[0])
            seizures_df = pd.DataFrame({'SEIZURES': seizures_values}, index=seizures_index)

            # Metas of Seizures in Dataframe
            meta_df = pd.read_csv(include_meta, delimiter=';')
            meta_df = meta_df[meta_df['UKLEEG_NUMBER'].isin([self.patient_id])]
            meta_df = meta_df[meta_df['TAG'].isin(seizures_values)].reset_index(drop=True)
            meta_df = meta_df.set_index(seizures_df.index)
            self.meta = meta_df.values.tolist()

            # Alle Spalten Droppen die für weitere Verarbeitung irrelavant sind
            meta_df = meta_df.drop(columns=['MYOCLONIC','EPILEPTIC','TYP I','TYP II','TYP III','SEITE','TAG','STUDY_ID','ORIGIN_START','EKG_QUAL','EEG','EMG','DELRE','DELLI','QUADRE','QUADLI','VIDEO','VIGILANZ','PERCEPTION','IMPULSE'])
            
            #print(self.feature_df.shape)
            self.feature_df = pd.concat([signal_df, seizures_df, meta_df, self.feature_df], ignore_index=False, axis=1)
            #print(self.feature_df.shape)
            #print(self.feature_df.columns)

        # Export File as Pickle
        if export is not None:
            self.feature_df.to_pickle(export)

    def plot_ecg_raw(self, sec_start=30, sec_pre=0, sec_post=20, save_graphic=None):
        '''
        Plots 10 seconds of a raw ecg signal 
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Data
        plt.plot(self.samples, label='Einthoven Lead II', linewidth=1.5, color='black')
        
        # Plot Bereich
        lower_limit = sec_start*self.sample_rate - sec_pre*self.sample_rate
        upper_limit = sec_start*self.sample_rate + sec_post * self.sample_rate
        
        # Plot Settings
        #plt.rcParams["font.family"] = "Times New Roman"
        #plt.rcParams["font.size"] = "12"
        plt.xlabel(f'Zeit[s]')
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-1, 3) 
        plt.ylabel('EKG [mV]')
        plt.grid(b=True, which='major', axis='both')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4) #fontsize='x-small'
        new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(-1 * sec_pre, sec_post, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.draw()
        plt.savefig(save_graphic+'EKG.png', dpi=300, format='png', transparent=False, bbox_inches='tight')
        
    def plot_ecg_preprocessed(self, sec_start=30, sec_pre=0, sec_post=20, save_graphic=None):
        '''
        Plots the different steps of the hamilton preprocessing for the 10 sec ecg signal 
        '''
        # Figur mit Subplots erstellen
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(12, 6)) 
        gs1 = gridspec.GridSpec(5, 1)
        gs1.update(hspace=0)

        # Plot Bereich
        lower_limit = sec_start*self.sample_rate - sec_pre*self.sample_rate
        upper_limit = sec_start*self.sample_rate + sec_post * self.sample_rate

        # Plot Data
        ax1 = plt.subplot(gs1[0])
        plt.axis('on')
        plt.plot(self.samples_detrended, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-1, 3)
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('Detrended')
      
        ax2 = plt.subplot(gs1[1])
        ax2.autoscale(enable=True, tight=True)
        plt.axis('on')
        plt.plot(self.samples_filtered, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('BP Filtered')
        
        ax3 = plt.subplot(gs1[2])
        ax3.autoscale(enable=True, tight=True)
        plt.axis('on')
        samples_diff = np.insert(self.samples_diff,0,self.samples_diff[0])
        plt.plot(samples_diff, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('Differentiated')
        
        ax4 = plt.subplot(gs1[3])
        ax4.autoscale(enable=True, tight=True)
        plt.axis('on')
        samples_rect = np.insert(self.samples_rect,0,self.samples_rect[0])
        plt.plot(samples_rect, label='Redctified', linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('Redctified')
        
        ax5 = plt.subplot(gs1[4])
        ax5.autoscale(enable=True, tight=True)
        plt.axis('on')
        preprocessed = np.insert(self.preprocessed,0,self.preprocessed[0])
        plt.plot(preprocessed, linewidth=1.5, color='black')
        plt.yticks([])
        plt.xlim(lower_limit, upper_limit)
        new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(-1 * sec_pre, sec_post, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.ylabel('MA Filtered')
        plt.xlabel(f'Zeit[s]')
        
        plt.draw()
        plt.savefig(save_graphic+'EKG_Preprocessed.png', dpi=300, format='png', transparent=False, bbox_inches='tight')
                  
    def plot_ecg_with_rpeaks(self, sec_start=30, sec_pre=0, sec_post=20, save_graphic=None):
        '''
        Plots 10 seconds of a raw ecg signal 
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Data
        plt.plot(self.samples, label='Einthoven Lead II', linewidth=1.5, color='black')
        peaks = [int(r_peak/self.period_ms) for r_peak in self.r_peaks]
        plt.plot(peaks, self.samples[peaks],'x',color='red', linewidth=4, label='Beats')

        # Plot Bereich
        lower_limit = sec_start*self.sample_rate - sec_pre*self.sample_rate
        upper_limit = sec_start*self.sample_rate + sec_post * self.sample_rate
        
        # Plot Settings
        plt.xlabel(f'Zeit[s]')
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-1, 3) 
        plt.ylabel('EKG [mV]')
        plt.grid(b=True, which='major', axis='both')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4) #fontsize='x-small'
        new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(-1 * sec_pre, sec_post, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.draw()
        plt.savefig(save_graphic + 'EKG_Beats.png', dpi=300, format='png', transparent=False, bbox_inches='tight')
    
    def plot_rr_intervals(self, sec_start=30, sec_pre=0, sec_post=20, save_graphic=None):
        '''
        Plots 10 seconds of a raw ecg signal 
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Bereich
        #lower_limit = sec_start * self.sample_rate - sec_pre * self.sample_rate
        lower_limit = (sec_start- sec_pre)*1000
        #upper_limit = sec_start * self.sample_rate + sec_post * self.sample_rate
        upper_limit = (sec_start+ sec_post)*1000
        
        # Grab Data
        beats = [beat for beat in self.r_peaks if beat in range(lower_limit, upper_limit)]
        ibi_interval = []
        for k in range(1,len(beats)):
            ibi_interval.append(math.ceil(beats[k] - beats[k - 1]))
        
        # Plot Data
        plt.bar(np.arange(len(ibi_interval)),ibi_interval, width=1, align='center',label='IBI Bars', color='white',edgecolor='black',linewidth=2)
        plt.plot(ibi_interval, 'x', color='red', linewidth=1.5)
        plt.plot(ibi_interval, color='red', linewidth=1.5, label='IBI Line')
        #plt.plot(peaks, self.samples[peaks], 'x', color='red', linewidth=4, label='Beats')
        
        # Plot Settings
        plt.xlabel(f'Intervall [k]')
        #plt.xlim(lower_limit, upper_limit)
        #plt.ylim(-1, 3) 
        plt.ylabel('Intervalllänge [ms]')
        #plt.grid(b=True, which='major', axis='both')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4)#fontsize='x-small'
        #new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        #plt.gca().set_xticks(new_xticks)
        new_time = [0,1,2,3,4,5,6,7,8]
        plt.gca().set_xticklabels(new_time)
        plt.draw()
        plt.savefig(save_graphic + 'RR_Intervals.png', dpi=300, format='png', transparent=False, bbox_inches='tight')
    
    def plot_rr_tachogram(self, sec_start=30, sec_pre=0, sec_post=20, save_graphic=None):
        '''
        Plots 10 seconds of a raw ecg signal 
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Bereich
        #lower_limit = sec_start * self.sample_rate - sec_pre * self.sample_rate
        lower_limit = (sec_start- sec_pre)*1000
        #upper_limit = sec_start * self.sample_rate + sec_post * self.sample_rate
        upper_limit = (sec_start+ sec_post)*1000
        
        # Grab Data
        beats = [beat for beat in self.r_peaks if beat in range(lower_limit, upper_limit)]
        ibi_interval = []
        for k in range(1,len(beats)):
            ibi_interval.append(math.ceil(beats[k] - beats[k - 1]))
        
        #n, bins = np.histogram(ibi_interval, bins=10)
        bin_sequence = np.arange(800,1400+1,50)
        plt.hist(ibi_interval,bins=bin_sequence,label='IBI',color='white',edgecolor='black',linewidth=2)
        # Plot Data
        # plt.bar(np.arange(len(ibi_interval)),ibi_interval, width=1, align='center',label='IBI Bars', color='white',edgecolor='black')
        # plt.plot(ibi_interval, 'x', color='red', linewidth=1.5)
        # plt.plot(ibi_interval, color='red', linewidth=1.5, label='IBI Line')
        #plt.plot(peaks, self.samples[peaks], 'x', color='red', linewidth=4, label='Beats')
        
        # Plot Settings
        plt.xlabel(f'Intervalllänge [ms]')
        #plt.xlim(lower_limit, upper_limit)
        #plt.ylim(-1, 3) 
        plt.ylabel('Absolute Häufigkeit')
        #plt.grid(b=True, which='major', axis='both')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4)#fontsize='x-small'
        #new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        #plt.gca().set_xticks(new_xticks)
        # new_time = [0,1,2,3,4,5,6,7,8]
        # plt.gca().set_xticklabels(new_time)
        plt.draw()
        plt.savefig(save_graphic + 'Tachogram.png', dpi=300, format='png', transparent=False, bbox_inches='tight')
     
        


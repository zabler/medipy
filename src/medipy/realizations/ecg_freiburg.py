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
from astropy.timeseries import LombScargle

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
        annotations = header['annotations']  #format Seconds.Milliseconds,Duration in sec,Tag: SEIZURE_P##_S##

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
        
        # Patient ID
        filename = path.split('/')[-1]
        self.patient_id = filename.split('_')[0]
            
    def data_cutter(self, duration_s=None):
        if duration_s is not None:
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
            
            # Merge all Dataframes
            self.feature_df = pd.concat([signal_df, seizures_df, meta_df, self.feature_df], ignore_index=False, axis=1)

        # Export File as Pickle
        if export is not None:
            self.feature_df.to_pickle(export)

    def plot_ecg_raw(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plot 10 seconds of a raw ecg signal 
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)

        # Plot Data
        plt.plot(self.samples, label='Einthoven Lead II', linewidth=1.5, color='black')
        
        # Plot Bereich
        lower_limit = start_sec_abs*self.sample_rate
        upper_limit = start_sec_abs*self.sample_rate + duration_sec_rel * self.sample_rate
        
        # Plot Settings
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-1, 3)
        plt.xlabel(f'Zeit[s]')
        plt.ylabel('EKG [mV]')
        plt.grid(b=True, which='major', axis='both')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4)
        new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic+'5210_Beispiel_Rohes_EKG_Signal.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
    def plot_ecg_preprocessed(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plots the preprocessing steps of a 10 seconds signal
        '''
        # Figur mit Subplots erstellen
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(12, 6))
        gs1 = gridspec.GridSpec(5, 1)
        gs1.update(hspace=0)

        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)

        # Plot Bereich
        lower_limit = start_sec_abs*self.sample_rate
        upper_limit = start_sec_abs*self.sample_rate + duration_sec_rel * self.sample_rate

        # Plot Data
        ax1 = plt.subplot(gs1[0])
        ax1.grid(True)
        ax1.set_xticklabels([])
        plt.plot(self.samples, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-1, 3)
        plt.yticks([])
        plt.ylabel('Raw Signal')

        ax2 = plt.subplot(gs1[1])
        ax2.grid(True)
        ax2.set_xticklabels([])
        plt.plot(self.samples_filtered, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.ylabel('BP Filtered')
        
        ax3 = plt.subplot(gs1[2])
        ax3.grid(True)
        ax3.set_xticklabels([])
        samples_diff = np.insert(self.samples_diff,0,self.samples_diff[0])
        plt.plot(samples_diff, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.ylabel('Differentiated')
        
        ax4 = plt.subplot(gs1[3])
        ax4.grid(True)  
        ax4.set_xticklabels([])
        samples_rect = np.insert(self.samples_rect, 0, self.samples_rect[0])
        plt.plot(samples_rect, label='Rectified', linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.ylabel('Rectified')
        
        ax5 = plt.subplot(gs1[4])
        ax5.grid(True)
        preprocessed = np.insert(self.preprocessed, 0, self.preprocessed[0])
        plt.plot(preprocessed, linewidth=1.5, color='black')
        plt.yticks([])
        plt.ylabel('MA Filtered')
        
        # Plot Settings
        plt.xlim(lower_limit, upper_limit)
        new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xlabel(f'Zeit[s]')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic+'5221_EKG_Preprocessing_des_Beispielsignals_nach_Hamilton.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
                  
    def plot_ecg_detected_rpeaks(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plots 10 seconds of a raw ecg signal and detected r peaks
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)

        # Plot Data
        plt.plot(self.samples, label='Einthoven Lead II', linewidth=1.5, color='black')
        peaks = [int(r_peak/self.period_ms) for r_peak in self.r_peaks]
        plt.plot(peaks, self.samples[peaks],'x',color=red, linewidth=4, label='R-Zacken')

        # Plot Bereich
        lower_limit = start_sec_abs*self.sample_rate
        upper_limit = start_sec_abs*self.sample_rate + duration_sec_rel * self.sample_rate
        
        # Plot Settings
        plt.xlabel(f'Zeit[s]')
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-1, 3) 
        plt.ylabel('EKG [mV]')
        plt.grid(b=True, which='major', axis='both')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4) 
        new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5222_Beispielsignal_mit_detektierten_R-Zacken.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
    
    def plot_rr_interval_tachogram(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        Plots the tachogram of detected in 10 seconds ecg signal
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)

        # Plot Bereich
        lower_limit = start_sec_abs*1000
        upper_limit = start_sec_abs*1000+duration_sec_rel*1000
        
        # Grab Data
        local_r_peaks = [r_peak for r_peak in self.r_peaks if r_peak in range(lower_limit, upper_limit)]
        local_rr_intervals = []
        for k in range(1, len(local_r_peaks)):
            local_rr_intervals.append(math.ceil(local_r_peaks[k] - local_r_peaks[k - 1]))
        
        # Plot Data
        plt.bar(np.arange(len(local_rr_intervals)), local_rr_intervals, width=1, align='center', label='RR-Intervalle', color='white', edgecolor='black', linewidth=1.5)
        plt.plot(local_rr_intervals, color=red, linewidth=1.5, label='RR-Intervall-Liste')
        
        # Plot Settings
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc='lower left', mode='expand', borderaxespad=0, ncol=4)
        plt.xlabel(f'Intervall [k]')
        plt.ylabel('Intervalllänge [ms]')
        new_time = [0, 1, 2, 3, 4, 5, 6, 7, 8] # Wegen 9 R-Zacken
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0,1400)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Rythmogramm_des_Beispielsignals.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
    
    def plot_rr_interval_histogram(self, start_sec_abs=30, duration_sec_rel=300, save_graphic=None):
        '''
        plots the histogram of rr_intervals of detected r peaks of 5min signal
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)

        # Plot Bereich
        lower_limit = start_sec_abs*1000
        upper_limit = start_sec_abs*1000+duration_sec_rel*1000
        
        # Grab Data
        local_r_peaks = [r_peak for r_peak in self.r_peaks if r_peak in range(lower_limit, upper_limit)]
        local_rr_intervals = []
        for k in range(1, len(local_r_peaks)):
            local_rr_intervals.append(math.ceil(local_r_peaks[k] - local_r_peaks[k - 1]))
        
        # Plot Data
        bin_sequence = np.arange(800,1400+1,25)
        plt.hist(local_rr_intervals, bins=bin_sequence, label='RR-Intervalle', color='white', edgecolor='black', linewidth=1.5)
               
        # Plot Settings
        plt.xlabel(f'Intervalllänge [ms]')
        plt.ylabel('Absolute Häufigkeit')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4)#fontsize='x-small'
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5224_Tachogramm_eines_5-Minuten-Abschnitts.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
     
    def plot_rr_interval_errors(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plots all rr intervall erros incl. the ecg origin for the exmaple signal
        '''
        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)

        # Plot Bereich
        lower_limit_s = start_sec_abs*self.sample_rate
        upper_limit_s = start_sec_abs*self.sample_rate + duration_sec_rel * self.sample_rate
        lower_limit = start_sec_abs*1000
        upper_limit = start_sec_abs*1000+duration_sec_rel*1000
        
        # Grab Origin Data
        local_samples = self.samples[lower_limit_s:upper_limit_s]
        local_r_peaks = [r_peak for r_peak in self.r_peaks if r_peak in range(lower_limit, upper_limit)]
        
        # ERROR TYPE A1: ECTOPIC LONG SHORT
        error_peaks = list(local_r_peaks)
        error_peaks[4] = error_peaks[4]+600
        error_rr_intervals = []
        for k in range(1, len(error_peaks)):
            error_rr_intervals.append(math.ceil(error_peaks[k] - error_peaks[k - 1]))

        # Plot Error Origin
        fig = plt.figure(figsize=(9, 2.5))
        plt.plot(local_samples, linewidth=1.5, color='black')
        ecg_error_peaks = [int((error_peak-lower_limit)/self.period_ms) for error_peak in error_peaks]
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks],'x',color=red, linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=red, linewidth=2)
        plt.ylim(-1, 3)
        plt.xlim(0, upper_limit_s - lower_limit_s)
        plt.grid(b=True, which='major', axis='both')
        new_xticks = np.arange(0, upper_limit_s-lower_limit_s + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Origin_A1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3,4], error_rr_intervals[3:5], width=1, align='center', color='white', edgecolor=wine, linewidth=3)
        plt.plot(error_rr_intervals, color=red, linewidth=2)
        new_xticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticks(new_xticks)
        new_time = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0, 2000)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Type_A1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # ERROR TYPE A2: ECTOPIC SHORT LONG
        error_peaks = list(local_r_peaks)
        error_peaks[4] = error_peaks[4]-600
        error_rr_intervals = []
        for k in range(1, len(error_peaks)):
            error_rr_intervals.append(math.ceil(error_peaks[k] - error_peaks[k - 1]))

        # Plot Error Origin
        fig = plt.figure(figsize=(9, 2.5))
        plt.plot(local_samples, linewidth=1.5, color='black')
        ecg_error_peaks = [int((error_peak-lower_limit)/self.period_ms) for error_peak in error_peaks]
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks],'x',color=red, linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=red, linewidth=2)
        plt.ylim(-1, 3)
        plt.xlim(0, upper_limit_s - lower_limit_s)
        plt.grid(b=True, which='major', axis='both')
        new_xticks = np.arange(0, upper_limit_s-lower_limit_s + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Origin_A2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3,4], error_rr_intervals[3:5], width=1, align='center', color='white', edgecolor=wine, linewidth=3)
        plt.plot(error_rr_intervals, color=red, linewidth=2)
        new_xticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticks(new_xticks)
        new_time = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0, 2000)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Type_A2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # ERROR TYPE B1: LONG INTERVALL
        error_peaks = list(local_r_peaks)
        error_peaks[4] = error_peaks[4]+254
        error_rr_intervals = []
        for k in range(1, len(error_peaks)):
            error_rr_intervals.append(math.ceil(error_peaks[k] - error_peaks[k - 1]))

        # Plot Error Origin
        fig = plt.figure(figsize=(9, 2.5))
        plt.plot(local_samples, linewidth=1.5, color='black')
        ecg_error_peaks = [int((error_peak-lower_limit)/self.period_ms) for error_peak in error_peaks]
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks],'x',color=red, linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=red, linewidth=2)
        plt.ylim(-1, 3)
        plt.xlim(0, upper_limit_s - lower_limit_s)
        plt.grid(b=True, which='major', axis='both')
        new_xticks = np.arange(0, upper_limit_s-lower_limit_s + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Origin_B1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3], error_rr_intervals[3], width=1, align='center', color='white', edgecolor=wine, linewidth=3)
        plt.plot(error_rr_intervals, color=red, linewidth=2)
        new_xticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticks(new_xticks)
        new_time = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0, 2000)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Type_B1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # ERROR TYPE B2: SHORT INTERVALL
        error_peaks = list(local_r_peaks)
        error_peaks[3] = error_peaks[3]-170
        error_rr_intervals = []
        for k in range(1, len(error_peaks)):
            error_rr_intervals.append(math.ceil(error_peaks[k] - error_peaks[k - 1]))

        # Plot Error Origin
        fig = plt.figure(figsize=(9, 2.5))
        plt.plot(local_samples, linewidth=1.5, color='black')
        ecg_error_peaks = [int((error_peak-lower_limit)/self.period_ms) for error_peak in error_peaks]
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks],'x',color=red, linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=red, linewidth=2)
        plt.ylim(-1, 3)
        plt.xlim(0, upper_limit_s - lower_limit_s)
        plt.grid(b=True, which='major', axis='both')
        new_xticks = np.arange(0, upper_limit_s-lower_limit_s + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Origin_B2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([2], error_rr_intervals[2], width=1, align='center', color='white', edgecolor=wine, linewidth=3)
        plt.plot(error_rr_intervals, color=red, linewidth=2)
        new_xticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticks(new_xticks)
        new_time = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0, 2000)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Type_B2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # ERROR TYPE C: EXTRA
        error_peaks = list(local_r_peaks)
        error_peaks.insert(5,error_peaks[4]+650)
        error_rr_intervals = []
        for k in range(1, len(error_peaks)):
            error_rr_intervals.append(math.ceil(error_peaks[k] - error_peaks[k - 1]))

        # Plot Error Origin
        fig = plt.figure(figsize=(9, 2.5))
        plt.plot(local_samples, linewidth=1.5, color='black')
        ecg_error_peaks = [int((error_peak-lower_limit)/self.period_ms) for error_peak in error_peaks]
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks],'x',color=red, linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=red, linewidth=2)
        plt.ylim(-1, 3)
        plt.xlim(0, upper_limit_s - lower_limit_s)
        plt.grid(b=True, which='major', axis='both')
        new_xticks = np.arange(0, upper_limit_s-lower_limit_s + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Origin_C.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([4,5], error_rr_intervals[4:6], width=1, align='center', color='white', edgecolor=wine, linewidth=3)
        plt.plot(error_rr_intervals, color=red, linewidth=2)
        new_xticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticks(new_xticks)
        new_time = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0, 2000)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Type_C.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # ERROR TYPE D: MISSED
        error_peaks = list(local_r_peaks)
        error_peaks = error_peaks[:4]+error_peaks[6:]
        error_rr_intervals = []
        for k in range(1, len(error_peaks)):
            error_rr_intervals.append(math.ceil(error_peaks[k] - error_peaks[k - 1]))

        # Plot Error Origin
        fig = plt.figure(figsize=(9, 2.5))
        plt.plot(local_samples, linewidth=1.5, color='black')
        ecg_error_peaks = [int((error_peak-lower_limit)/self.period_ms) for error_peak in error_peaks]
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks],'x',color=red, linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=red, linewidth=2)
        plt.ylim(-1, 3)
        plt.xlim(0, upper_limit_s - lower_limit_s)
        plt.grid(b=True, which='major', axis='both')
        new_xticks = np.arange(0, upper_limit_s-lower_limit_s + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Origin_D.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        
        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3], error_rr_intervals[3], width=1, align='center', color='white', edgecolor=wine, linewidth=3)
        plt.plot(error_rr_intervals, color=red, linewidth=2)
        new_xticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticks(new_xticks)
        new_time = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0, 2000)
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Error_Type_D.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def plot_psd(self, start_sec_abs=30, duration_sec_rel=300, save_graphic=None):
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Colors
        blue = (0, 0.4470, 0.7410)
        red = (0.8500, 0.3250, 0.0980)
        yellow = (0.9290, 0.6940, 0.1250)
        purple = (0.4940, 0.1840, 0.5560)
        grey = (0.5140, 0.5140, 0.5140)
        wine = (0.6350, 0.0780, 0.1840)

        # Plot Bereich
        lower_limit = start_sec_abs*1000
        upper_limit = start_sec_abs*1000+duration_sec_rel*1000
        
        # Grab Data
        rri_sum = 0
        rr_intervals = []
        for rri in self.rr_intervals:
            rri_sum = rri_sum + rri
            if rri_sum in range(lower_limit, upper_limit):
                rr_intervals.append(rri)
        #rr_intervals = [rri for rri in self.rr_intervals if rri in range(lower_limit, upper_limit)]

        # Data
        nni_tmstp = np.cumsum(rr_intervals) / 100
        timestamp_list = nni_tmstp - nni_tmstp[0]

        freq, psd = LombScargle(timestamp_list, rr_intervals,
                                normalization='psd').autopower(minimum_frequency=0.04,
                                                               maximum_frequency=0.40)
        
        # Plot
        plt.plot(freq, psd)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5224_PSD.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')
        


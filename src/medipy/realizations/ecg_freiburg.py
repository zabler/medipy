'''
ecg_freiburg.py
'''

import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from pyedflib import highlevel
from astropy.timeseries import LombScargle
import nolds
import csv  
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
        self.colours = {'blue': (0, 0.4470, 0.7410), 'red': (0.8500, 0.3250, 0.0980), 'yellow': (0.9290, 0.6940, 0.1250), 'purple': (0.4940, 0.1840, 0.5560), 'grey': (0.5140, 0.5140, 0.5140), 'wine': (0.6350, 0.0780, 0.1840)}

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
        '''
        This method cuts the samples of the signals and its accordings values to the given duration
        '''
        if duration_s is not None:
            stop = duration_s * self.sample_rate
            self.samples = np.array(self.samples[0:stop])
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
        time_feature_names = ['RR_MEDIAN', 'RR_MEAN', 'SDNN', 'SDSD', 'RMSSD', 'NN50', 'pNN50', 'HR_MEAN', 'HR_MAX_MIN']
        frequency_feature_names = ['LF_POWER', 'HF_POWER', 'LF_HF_RATIO', 'LF_NU', 'HF_NU', 'LF_PEAK', 'HF_PEAK']
        nonlinear_feature_names = ['SD1', 'SD2', 'CSI', 'MODIFIED_CSI', 'CVI', 'DF_ALPHA_1', 'DF_ALPHA_2']

        # Samples in Dataframe
        samples_df = pd.DataFrame({'SAMPLES': self.samples}, index=self.grid)

        # R-Peaks, RR_Intervals and RR_Missings to Dataframe
        rr_intervals = self.rr_intervals.copy()
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
            if self.rr_plausibility_check(rr_intervals_window, window=300, normal_level=0.1, artefact_level=0.01):
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
            meta_df = meta_df.drop(columns=['MYOCLONIC', 'EPILEPTIC', 'TYP I', 'TYP II', 'TYP III', 'SEITE', 'TAG', 'STUDY_ID', 'ORIGIN_START', 'EKG_QUAL', 'EEG', 'EMG', 'DELRE', 'DELLI', 'QUADRE', 'QUADLI', 'VIDEO', 'VIGILANZ', 'PERCEPTION', 'IMPULSE'])

            # Merge all Dataframes
            self.feature_df = pd.concat([signal_df, seizures_df, meta_df, self.feature_df], ignore_index=False, axis=1)

        # Export File as Pickle
        if export is not None:
            self.feature_df.to_pickle(export)

    def short_term_hrv_marked_extractor(self, window=300, overlap=0.999, area_min=24, include_meta=None, export=None):
        '''
        Short Term HRV Marked Extraktor
        Iteration mittels geleitendem (Overlap=1) Fenster / oder springendem (Overlap<1) Fenster über Grid (entspricht Index des Dataframes),
        wenn ein Anfall in einer Area+- markiert ist
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
        time_feature_names = ['RR_MEDIAN', 'RR_MEAN', 'SDNN', 'SDSD', 'RMSSD', 'NN50', 'pNN50', 'HR_MEAN', 'HR_MAX_MIN']
        frequency_feature_names = ['LF_POWER', 'HF_POWER', 'LF_HF_RATIO', 'LF_NU', 'HF_NU', 'LF_PEAK', 'HF_PEAK']
        nonlinear_feature_names = ['SD1', 'SD2', 'CSI', 'MODIFIED_CSI', 'CVI', 'DF_ALPHA_1', 'DF_ALPHA_2']

        # Samples in Dataframe
        samples_df = pd.DataFrame({'SAMPLES': self.samples}, index=self.grid)

        # R-Peaks, RR_Intervals and RR_Missings to Dataframe
        rr_intervals = self.rr_intervals.copy()
        rr_intervals = np.insert(rr_intervals, 0, int(np.mean(self.rr_intervals[0:10])))
        r_peaks_df = pd.DataFrame({'R_PEAKS': np.array(self.r_peaks), 'RR_INTERVALS': rr_intervals}, index=np.array(self.r_peaks))

        # Combine in Signal Dataframe
        signal_df = pd.concat([samples_df, r_peaks_df], ignore_index=False, axis=1).fillna(value=np.nan)
        #signal_df = signal_df.drop(columns=['SAMPLES', 'RR_INTERVALS'])
        print(signal_df.shape)

        # Get from Signal Dataframe including time information
        rr_intervals_frame = signal_df['RR_INTERVALS'].to_numpy()

        # Seizures in Dataframe
        seizures_values = []
        seizures_index = []
        for seizure in self.seizures:
            seizures_values.append(seizure[1])
            seizures_index.append(seizure[0])
        seizures_df = pd.DataFrame({'SEIZURES': seizures_values}, index=seizures_index)
        index_check = seizures_df.index.to_numpy()
        
        # Create Feature Dataframe
        self.feature_df = pd.DataFrame(np.nan, index=self.grid, columns=time_feature_names + frequency_feature_names + nonlinear_feature_names)
        print(self.feature_df.shape)
    
        # HRV Calculation # BEDINGUNG SEIZRUE IN WINDOW AREA
        half_window = int((window / 2) * self.sample_rate)
        step_size = int(round((1 - overlap) * half_window * 2))
        area = int(area_min/2)*60*self.sample_rate
        if step_size == 0:
            step_size = 1
        steps = np.arange(0, len(self.grid), step=step_size)
        for step in steps:
            if step < half_window or step+half_window > steps[-1]:
                continue
            area_checks = abs(index_check-step*self.period_ms)
            if not any(area_check < area for area_check in area_checks):
                continue
            rr_intervals_window = rr_intervals_frame[step - half_window:step + half_window]
            rr_intervals_window = rr_intervals_window[~np.isnan(rr_intervals_window)]
            if self.rr_plausibility_check(rr_intervals_window, window=300, normal_level=0.1, artefact_level=0.01):
                time_features = self.hrv_features_time(rr_intervals_window)
                frequency_features = self.hrv_features_frequency(rr_intervals_window)
                nonlinear_features = self.hrv_features_nonlinear(rr_intervals_window)
                features = {**time_features, **frequency_features, **nonlinear_features}
                for feature in features:
                    self.feature_df.at[self.grid[step], feature] = features[feature]
            
            # andere Option: Kein Durchlauf des Gesamtsignals, nur um Seizure herum, Tabelle(n) werden für Seizures gespeichert, Analyse muss angepasst werden
            # weiter: Signal_df nicht abspeichern oder signal_df abspeichern aber Signal und R Peaks droppen? oder doch dazu speichern

        # Append Meta Infos
        if include_meta is not None:

            # Metas of Seizures in Dataframe
            meta_df = pd.read_csv(include_meta, delimiter=';')
            meta_df = meta_df[meta_df['UKLEEG_NUMBER'].isin([self.patient_id])]
            meta_df = meta_df[meta_df['TAG'].isin(seizures_values)].reset_index(drop=True)
            meta_df = meta_df.set_index(seizures_df.index)
            self.meta = meta_df.values.tolist()

            # Alle Spalten Droppen die für weitere Verarbeitung irrelavant sind
            meta_df = meta_df.drop(columns=['MYOCLONIC', 'EPILEPTIC', 'TYP I', 'TYP II', 'TYP III', 'SEITE', 'TAG', 'STUDY_ID', 'ORIGIN_START', 'EKG_QUAL', 'EEG', 'EMG', 'DELRE', 'DELLI', 'QUADRE', 'QUADLI', 'VIDEO', 'VIGILANZ', 'PERCEPTION', 'IMPULSE'])

            # Merge
            self.feature_df = pd.concat([signal_df, seizures_df, meta_df, self.feature_df], ignore_index=False, axis=1) #
        else:
            # Merge without Meta
            self.feature_df = pd.concat([signal_df, seizures_df, self.feature_df], ignore_index=False, axis=1) #

        # Export File as Pickle
        if export is not None:
            self.feature_df.to_pickle(export)

    def plot_ecg_raw(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plots a raw ecg signal
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Data
        plt.plot(self.samples, label='Einthoven Lead II', linewidth=1.5, color='black')

        # Plot Bereich
        lower_limit = start_sec_abs*self.sample_rate
        upper_limit = start_sec_abs*self.sample_rate + duration_sec_rel * self.sample_rate

        # Plot Settings
        plt.xlim(lower_limit, upper_limit)
        plt.ylim(-1, 3)
        plt.xlabel('Zeit[s]')
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
        plots preprocessing steps
        '''
        # Figur mit Subplots erstellen
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(12, 6))
        gs1 = gridspec.GridSpec(5, 1)
        gs1.update(hspace=0)

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
        plt.ylabel('Rohsignal')

        ax2 = plt.subplot(gs1[1])
        ax2.grid(True)
        ax2.set_xticklabels([])
        plt.plot(self.samples_filtered, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.ylabel('BP-Filter')

        ax3 = plt.subplot(gs1[2])
        ax3.grid(True)
        ax3.set_xticklabels([])
        samples_diff = np.insert(self.samples_diff, 0, self.samples_diff[0])
        plt.plot(samples_diff, linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.ylabel('5P-Diff')

        ax4 = plt.subplot(gs1[3])
        ax4.grid(True)
        ax4.set_xticklabels([])
        samples_rect = np.insert(self.samples_rect, 0, self.samples_rect[0])
        plt.plot(samples_rect, label='Rectified', linewidth=1.5, color='black')
        plt.xlim(lower_limit, upper_limit)
        plt.yticks([])
        plt.ylabel('Gleichrichter')

        ax5 = plt.subplot(gs1[4])
        ax5.grid(True)
        preprocessed = np.insert(self.preprocessed, 0, self.preprocessed[0])
        plt.plot(preprocessed, linewidth=1.5, color='black')
        plt.yticks([])
        plt.ylabel('GM-Filter')

        # Plot Settings
        plt.xlim(lower_limit, upper_limit)
        new_xticks = np.arange(lower_limit, upper_limit + 1, 2 * self.sample_rate)
        plt.gca().set_xticks(new_xticks)
        new_time = np.linspace(0, duration_sec_rel, len(new_xticks), endpoint=True, dtype=str)
        plt.gca().set_xticklabels(new_time)
        plt.xlabel('Zeit[s]')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic+'5221_EKG_Preprocessing_des_Beispielsignals_nach_Hamilton.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def plot_ecg_detected_rpeaks(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plots a raw ecg signal and detected r peaks
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Data
        plt.plot(self.samples, label='Einthoven Lead II', linewidth=1.5, color='black')
        peaks = [int(r_peak/self.period_ms) for r_peak in self.r_peaks]
        plt.plot(peaks, self.samples[peaks], 'x', color=self.colours['wine'], linewidth=4, label='R-Zacken', markersize=8)

        # Plot Bereich
        lower_limit = start_sec_abs*self.sample_rate
        upper_limit = start_sec_abs*self.sample_rate + duration_sec_rel * self.sample_rate

        # Plot Settings
        plt.xlabel('Zeit[s]')
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

    def plot_rr_interval_bartachogram(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plots the tachogram in bars and as a line of rr intervals
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

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
        plt.plot(local_rr_intervals, color=self.colours['blue'], linewidth=1.5, label='RR-Intervall-Folge')

        # Plot Settings
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc='lower left', mode='expand', borderaxespad=0, ncol=4)
        plt.xlabel('Intervall [k]')
        plt.ylabel('Intervalllänge [ms]')
        new_time = [0, 1, 2, 3, 4, 5, 6, 7, 8] # Wegen 9 R-Zacken
        plt.gca().set_xticklabels(new_time)
        plt.ylim(0, 1400)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5223_Tachogramm_des_Beispielsignals.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def plot_rr_interval_errors(self, start_sec_abs=30, duration_sec_rel=10, save_graphic=None):
        '''
        plots all rr intervall erros incl. the ecg origin for the example signal
        '''
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
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks], 'x', color=self.colours['wine'], linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=self.colours['wine'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Origin_A1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3, 4], error_rr_intervals[3:5], width=1, align='center', color='white', edgecolor=self.colours['wine'], linewidth=3)
        plt.plot(error_rr_intervals, color=self.colours['blue'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Type_A1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

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
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks], 'x', color=self.colours['wine'], linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=self.colours['wine'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Origin_A2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3, 4], error_rr_intervals[3:5], width=1, align='center', color='white', edgecolor=self.colours['wine'], linewidth=3)
        plt.plot(error_rr_intervals, color=self.colours['blue'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Type_A2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

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
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks], 'x', color=self.colours['wine'], linewidth=3, markersize=8, label='R-Zacken')
        #plt.vlines(ecg_error_peaks, ymin=-1,ymax=3,color=self.colours['wine'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Origin_B1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3], error_rr_intervals[3], width=1, align='center', color='white', edgecolor=self.colours['wine'], linewidth=3)
        plt.plot(error_rr_intervals, color=self.colours['blue'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Type_B1.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

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
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks], 'x', color=self.colours['wine'], linewidth=3, markersize=8, label='R-Zacken')
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
            plt.savefig(save_graphic + '5224_Error_Origin_B2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([2], error_rr_intervals[2], width=1, align='center', color='white', edgecolor=self.colours['wine'], linewidth=3)
        plt.plot(error_rr_intervals, color=self.colours['blue'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Type_B2.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # ERROR TYPE C: EXTRA
        error_peaks = list(local_r_peaks)
        error_peaks.insert(5, error_peaks[4] + 650)
        error_rr_intervals = []
        for k in range(1, len(error_peaks)):
            error_rr_intervals.append(math.ceil(error_peaks[k] - error_peaks[k - 1]))

        # Plot Error Origin
        fig = plt.figure(figsize=(9, 2.5))
        plt.plot(local_samples, linewidth=1.5, color='black')
        ecg_error_peaks = [int((error_peak-lower_limit)/self.period_ms) for error_peak in error_peaks]
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks], 'x', color=self.colours['wine'], linewidth=3, markersize=8, label='R-Zacken')
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
            plt.savefig(save_graphic + '5224_Error_Origin_C.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([4, 5], error_rr_intervals[4:6], width=1, align='center', color='white', edgecolor=self.colours['wine'], linewidth=3)
        plt.plot(error_rr_intervals, color=self.colours['blue'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Type_C.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

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
        plt.plot(ecg_error_peaks, local_samples[ecg_error_peaks], 'x', color=self.colours['wine'], linewidth=3, markersize=8, label='R-Zacken')
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
            plt.savefig(save_graphic + '5224_Error_Origin_D.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

        # Plot Error Type
        fig = plt.figure(figsize=(9, 2.5))
        plt.bar(np.arange(len(error_rr_intervals)), error_rr_intervals, width=1, align='center', color='white', edgecolor='black', linewidth=1.5)
        plt.bar([3], error_rr_intervals[3], width=1, align='center', color='white', edgecolor=self.colours['wine'], linewidth=3)
        plt.plot(error_rr_intervals, color=self.colours['blue'], linewidth=2)
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
            plt.savefig(save_graphic + '5224_Error_Type_D.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def plot_rr_interval_histogram(self, start_sec_abs=30, duration_sec_rel=300, save_graphic=None):
        '''
        plots the histogram of rr_intervals of detected r peaks
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()

        # Plot Bereich
        lower_limit = start_sec_abs*1000
        upper_limit = start_sec_abs*1000+duration_sec_rel*1000

        # Grab Data
        local_r_peaks = [r_peak for r_peak in self.r_peaks if r_peak in range(lower_limit, upper_limit)]
        local_rr_intervals = []
        for k in range(1, len(local_r_peaks)):
            local_rr_intervals.append(math.ceil(local_r_peaks[k] - local_r_peaks[k - 1]))
        local_rr_median = np.median(local_rr_intervals)
        local_rr_mean = np.mean(local_rr_intervals)
        real_number = len(local_rr_intervals)
        theo_number = int((300 * 1000) / local_rr_median)
        theo2_number = int((300 * 1000) / local_rr_mean)

        # Plot Data
        bin_sequence = np.arange(800, 1400 + 1, 10)
        ax1.hist(local_rr_intervals, bins=bin_sequence, label='RR-Intervalle', color='white', edgecolor='black', linewidth=1.5)
        #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr = '\n'.join((f'Theoretische Anzahl RR-Intervalle im Median: {theo_number}', f'Theoretische Anzahl RR-Intervalle im Mittel: {theo2_number}', f'Tatsächliche Anzahl RR-Intervalle: {real_number}'))
        ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=8, verticalalignment='top')#, bbox=props)
        ax1.set_xlim(800, 1400)
        ax2.set_xlim(800, 1400)
        ax2.set_xticks([])
        ax2.axvline(x=local_rr_median, label='Median', color=self.colours['wine'], linewidth=4)
        ax2.axvline(x=local_rr_mean, label='Mittelwert', color=self.colours['blue'], linewidth=4)


        # Plot Settings
        ax1.set_xlabel('Intervalllänge [ms]')
        ax1.set_ylabel('Absolute Häufigkeit')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc="lower left", mode='expand', borderaxespad=0, ncol=4)#fontsize='x-small'
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '5224_Histogramm_eines_5-Minuten-Abschnitts_und_Intervallschätzung.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def plot_rr_interval_tachogram(self, start_sec_abs=30, duration_sec_rel=300, save_graphic=None):
        '''
        plots the tachogram only as a line of rr intervals
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Bereich
        lower_limit = start_sec_abs*1000
        upper_limit = start_sec_abs*1000+duration_sec_rel*1000

        # Grab Data
        local_r_peaks = [r_peak for r_peak in self.r_peaks if r_peak in range(lower_limit, upper_limit)]
        local_rr_intervals = []
        for k in range(1, len(local_r_peaks)):
            local_rr_intervals.append(math.ceil(local_r_peaks[k] - local_r_peaks[k - 1]))

        # Plot Data
        plt.plot(local_rr_intervals, color='black', linewidth=1.5, label='NN-Intervall-Folge')

        # Plot Settings
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc='lower left', mode='expand', borderaxespad=0, ncol=4)
        plt.xlabel('Intervall [k]')
        plt.ylabel('Intervalllänge [ms]')
        plt.ylim(800, 1400)
        plt.xlim(0, len(local_rr_intervals)-1)
        new_time = [1, 50, 100, 150, 200, 250]
        plt.gca().set_xticklabels(new_time)
        # plt.axis('off')
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '2211_Kontinuierliche_Darstellung_einer_NN-Intervall-Folge.svg', dpi=300, format='svg', transparent=True, bbox_inches='tight')

    def plot_rr_interval_pointcare(self, start_sec_abs=30, duration_sec_rel=300, save_graphic=None):
        '''
        plots the pointcare diagram of rr intervals
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

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
        x = rr_intervals.copy()
        x.pop(-1)
        y = rr_intervals.copy()
        y.pop(0)
        drr = np.diff(rr_intervals)
        sd1 = np.sqrt(np.std(drr, ddof=1) ** 2 * 0.5)
        sd2 = np.sqrt(2 * np.std(rr_intervals, ddof=1)** 2 - 0.5 * np.std(drr, ddof=1)** 2)

        # Plot
        ax.scatter(x, y, s=5, label='Poincare', color='black')
        na = 4
        ax.arrow(np.mean(rr_intervals), np.mean(rr_intervals), (-sd1 + na) * np.cos(np.deg2rad(45)), (sd1 - na) * np.sin(np.deg2rad(45)), head_width=5, head_length=5, fc=self.colours['blue'], ec=self.colours['blue'], zorder=4, linewidth=1.5)
        ax.arrow(np.mean(rr_intervals), np.mean(rr_intervals), (sd2 - na) * np.cos(np.deg2rad(45)), (sd2 - na) * np.sin(np.deg2rad(45)), head_width=5, head_length=5, fc=self.colours['wine'], ec=self.colours['wine'], zorder=4, linewidth=1.5)
        ax.add_line(mpl.lines.Line2D((np.mean(rr_intervals) - sd1 * np.cos(np.deg2rad(45)) * na, np.mean(rr_intervals) + sd1 * np.cos(np.deg2rad(45)) * na), (np.mean(rr_intervals) + sd1 * np.sin(np.deg2rad(45)) * na, np.mean(rr_intervals) - sd1 * np.sin(np.deg2rad(45)) * na), c=self.colours['blue'], ls=':', alpha=0.6))
        ax.add_line(mpl.lines.Line2D((800, 1400), (800, 1400), c=self.colours['wine'], ls=':', alpha=0.6))
        ax.set_xlim(800, 1400)
        ax.set_ylim(800, 1400)
        ax.scatter(1401, 1401, label='SD1', color=self.colours['blue'])
        ax.scatter(1401, 1401, label='SD2', color=self.colours['wine'])
        ax.set_xlabel(r'$NN_j [ms]$')
        ax.set_ylabel(r'$NN_{j+1} [ms]$')
        ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc='lower left', mode='expand', borderaxespad=0, ncol=4)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '2213_Poincare_Diagramm.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def plot_rr_interval_dfa(self, start_sec_abs=30, duration_sec_rel=300, save_graphic=None):
        '''
        plots the dfa diagram of rr intervals
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(6, 6))

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

        short = range(4, 16+1)
        long = range(17, 64+1)
        alpha1, dfa_short = nolds.dfa(rr_intervals, short, debug_data=True, overlap=False)
        vals1 = dfa_short[0]
        flucts1 = dfa_short[1]
        poly1 = np.polyval(dfa_short[2], dfa_short[0])
        alpha2, dfa_long = nolds.dfa(rr_intervals, long, debug_data=True, overlap=False)
        vals2 = dfa_long[0]
        flucts2 = dfa_long[1]
        poly2 = np.polyval(dfa_long[2], dfa_long[0])

        # Plot
        plt.plot(vals1, flucts1, 'o', color='black', markersize=3, label='DFA')
        plt.plot(vals1, poly1, color=self.colours['blue'], label=r'$\alpha_{1}$')#: '+f'{alpha1.round(decimals=2)}')
        plt.plot(vals2, flucts2, 'o', markersize=3, color='black')
        plt.plot(vals2, poly2, color=self.colours['wine'], label=r'$\alpha_{2}$')#: ' + f'{alpha2.round(decimals=2)}')

        # plt.xlim(1.0, 3.5)
        # plt.ylim(3.5, 6.0)
        plt.xlabel('ln [n]')
        plt.ylabel('ln F[n]')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc='lower left', mode='expand', borderaxespad=0, ncol=4)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '2213_Trendbereinigte_Fluktuationsanalyse.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def plot_rr_interval_psd(self, start_sec_abs=30, duration_sec_rel=300, save_graphic=None):
        '''
        plots the psd of rr intervals
        '''
        # Figur Erstellen
        fig = plt.figure(figsize=(12, 4))

        # Plot Bereich
        lower_limit = start_sec_abs*1000
        upper_limit = start_sec_abs*1000+duration_sec_rel*1000

        # Grab Data
        rri_sum = 0
        rr_intervals = []
        for rri in self.rr_intervals:
            rri_sum = rri_sum + rri
            if rri_sum in range(lower_limit, upper_limit):
                rr_intervals.append(int(rri))

        rr_timestamps_cumsum = np.cumsum(rr_intervals) /1000 # in sec damit Hz
        rr_timestamps = (rr_timestamps_cumsum - rr_timestamps_cumsum[0])
        rr_intervals = np.divide(rr_intervals, 1000)

        freq, psd = LombScargle(rr_timestamps, rr_intervals, normalization='psd').autopower(minimum_frequency=0.040, maximum_frequency=0.400, samples_per_peak=5)
        # psd = psd * 100000 # damit Magnitude in ms^2/Hz
        lf = [f for f in freq if f <= 0.150]
        hf = [f for f in freq if f > 0.150]

        # Plot
        plt.plot(freq, psd, label='Periodogramm', color='black', linewidth=1.5)
        plt.fill_between(lf, psd[:len(lf)], color=self.colours['blue'], label=r'$LF_{power}$')
        plt.fill_between(hf, psd[len(lf):], color=self.colours['wine'], label=r'$HF_{power}$')
        plt.axvline(x=lf[np.argmax(psd[:len(lf)])], color=self.colours['blue'], label=r'$LF_{peak}$', linestyle='dashed')
        plt.axvline(x=hf[np.argmax(psd[len(lf):])], color=self.colours['wine'], label=r'$HF_{peak}$', linestyle='dashed')
        plt.xlim(0.040, 0.400)
        plt.ylim(0, 0.3)
        plt.xlabel(r'Frequenz $[Hz]$')
        plt.ylabel(r'Magnitude $[s^2/Hz]$')
        plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.5), loc='lower left', mode='expand', borderaxespad=0, ncol=5)
        plt.draw()
        if save_graphic is not None:
            plt.savefig(save_graphic + '2212_Schätzung_des_Leistungsdichtespektrums_der_NN_Intervall_Folge.svg', dpi=300, format='svg', transparent=False, bbox_inches='tight')

    def export_error_stats(self, path, name):
        name = [name]

        window_stats = []
        window_stats.append(self.plausible+self.unplausible_no_data+self.unplausible_not_normal+self.unplausible_artefacts)
        window_stats.append(self.plausible)
        window_stats.append(self.unplausible_no_data+self.unplausible_not_normal+self.unplausible_artefacts)
        window_stats.append(self.unplausible_no_data)
        window_stats.append(self.unplausible_not_normal)
        window_stats.append(self.unplausible_artefacts)
        
        rr_stats = []
        artefacts_list = self.ectopic_intervals + self.long_short_intervals + self.extra_intervals + self.missed_intervals
        if self.rr_checked:
            rr_stats.append(int(np.cumsum(self.rr_checked)[-1]))
        else:
            rr_stats.append(0)
        if artefacts_list:
            rr_stats.append(int(np.cumsum(artefacts_list)[-1]))
        else:
            rr_stats.append(0)
        if self.ectopic_intervals:
            rr_stats.append(int(np.cumsum(self.ectopic_intervals)[-1]))
        else:
            rr_stats.append(0)
        if self.long_short_intervals:
            rr_stats.append(int(np.cumsum(self.long_short_intervals)[-1]))
        else:
            rr_stats.append(0)
        if self.extra_intervals:
            rr_stats.append(int(np.cumsum(self.extra_intervals)[-1]))
        else:
            rr_stats.append(0)
        if self.missed_intervals:
            rr_stats.append(int(np.cumsum(self.missed_intervals)[-1]))
        else:
            rr_stats.append(0)
        
        error_stats = name+window_stats+rr_stats

        with open(path, 'a') as file:
            writer = csv.writer(file)
            writer.writerow(error_stats)
    
    


        
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

        # tags = ['SEIZURE']
        # for annotation in annotations:
        #     self.tags.append([annotation_to_grid(self, annotation[0]), annotation[2]])
        #     if any(tag in annotation[2] for tag in tags):
        #         self.seizures.append([annotation_to_grid(self, annotation[0]), annotation[2]])
                
        tags = ['SEIZURE']
        for annotation in annotations:
            if any(tag in annotation[2] for tag in tags):
                self.seizures.append([annotation_to_grid(self, annotation[0]), annotation[2]])
            else:
                self.tags.append([annotation_to_grid(self, annotation[0]), annotation[2]])

        # tags = ['SEIZURE']
        # for annotation in annotations:
        #     if any(tag in annotation[2] for tag in tags):
        #         self.tags.append([annotation_to_grid(self, annotation[0]), annotation[2]])
    
    def short_time_hrv_calculator(self, window=2, overlap=1):
        '''
        st_hrv_calculator renntü ber eine PAndas table
        Berechnet HRV Feature XY, basierend auf Zeitbasis als Index mit zentralem Speicherzeitpunkt, RRList im Zeitraum, und markierter Fehler in einem Zeitfenster von 5min überlappend ?
        je nach FEature und je nach markierten Fehler 
        in ecg freiburg
        '''
        # Feature Names
        time_features = ['SDNN','SDRR','SDANN','SDNN Index','RMSSD','NN50','pNN50','DeltaHRMaxHRMin','TRI','TINN']
        freq_features = ['VLF_Power','LF_Power','LF_rel','LF_Nu','LF_Peak','HF_Power','HF_rel','HF_Nu','HF_Peak','LF/HF']
        non_linear_features = ['S', 'SD1', 'SD2', 'SD1/SD2', 'ApEn', 'SampEn', 'DF_alpha1', 'DF_alpha2', 'D_2']
        
        # Samples in Dataframe
        samples_df = pd.DataFrame({'SAMPLES': self.samples})

        # R-Peaks to Dataframe
        rr_intervals = self.rr_intervals.insert(0, np.mean(self.rr_intervals[0:10]))
        rr_missings = self.rr_missing.insert(0, 0)
        r_peaks_df = pd.DataFrame({'R_PEAKS': rr_intervals, 'RR_MISSING': rr_missings}, index=self.r_peaks)

        # Iterate Dataframe
        work_df = pd.concat([samples_df, r_peaks_df], ignore_index=False, axis=1)

        # Feature Dataframe
        feature_df = pd.DataFrame(index=work_df.index,columns=time_features+freq_features+non_linear_features) für alle Features dataframe mit nans und index time

        for work_index, row in work_df.iterrows(): # notwenig über diesen zu itterieren oder über range sinnvoller?# evtl nur über den index # pd besser da nicht kompletter satensatz evtl gelesen wird?
            if work_index < self.samplerate * int(window/2):
                break
            current_rr_interval = work_df.iloc[work_index-150:work_index*150,1] #ignore np.nan # work index vermutlich andere werte # nur mit pandas einfach einsammeln der zracken im hießigen bereich möglich
            current_rr_missing = work_df.iloc[work_index-150:work_index*150,2]
            if self.missing_data_check(current_rr_interval, current_rr_missing, güte=0.1):
                # alle funktionen nacheindnaer aufrufen, rr werte und misssing geben, interpolieren, feature berchnen und in feature tabelle einfügen mittelsindex
                feature_df.set_value(workd_index, 'SDANN', ecg.hrv_feature_sdann(ecg.rr_interpolation(method=linear)))

            # elif für andere güte durchführbar, bsp für 20% fehlende werte, andere interpolations methodik
            else:
                pass
                # für alle nan schrieben, da nan schon gegebn gehe weiter
            #work_index hochzählen anhand overlapping
            #damit nicht für jeden Sample etwas gerechnet wird 
                
        
        # über index lösen
        # HRV Extraktion

    def export_feature_table_with_classes(self, feature_table,file, meta_path):
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
        
        # Seizures in Dataframe
        tags_values = []
        tags_index = []
        for tag in self.tags:
            tags_values.append(tag[1])
            tags_index.append(tag[0])
        tags_df = pd.DataFrame({'TAGS': tags_values}, index=tags_index)

        # Metas of Seizures in Dataframe
        meta_df = pd.read_csv(meta_path, delimiter=';')
        meta_df = meta_df[meta_df['TAG'].isin(tags_values)].reset_index(drop=True)
        meta_df = meta_df.set_index(tags_df.index)
        self.meta = meta_df.values.tolist()

        # Alle Spalten Droppen die für weitere Verarbeitung irrelavant sind
        meta_df = meta_df.drop(columns=['MYOCLONIC','EPILEPTIC','TYP I','TYP II','TYP III','SEITE','TAG','STUDY_ID','ORIGIN_START','EMG','DELRE','DELLI','QUADRE','QUADLI','VIDEO','VIGILANZ','PERCEPTION','IMPULSE'])  
        
        # # Samples in Dataframe
        # samples_df = pd.DataFrame({'SAMPLES': self.samples})

        # # R-Peaks in Dataframe
        # rr_intervals = self.rr_intervals.insert(0, np.nan)
        # r_peaks_df = pd.DataFrame({'R_PEAKS': rr_intervals}, index=self.r_peaks)

        # Merge Dataframe
        all_df = pd.concat([feature_table, tags_df, meta_df], ignore_index=False, axis=1)

        #Save Dataframe to pickle (evtl. zu Testzwecken als csv)
        all_df.to_pickle(file)


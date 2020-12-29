'''
This is the hr_features.py file
'''
import pandas as pd
from medipy.interfaces.features import Features

class Hrv(pd.DataFrame, Features): #Erben von der Klasse und oder der Instanz, Ist Feautres automatisch mit implementiert??
    '''
    This is the hr_features class
    '''
    # def __init__(self):
    #     super(HrFeatures, self).__init__()

    # def create_index(self, length, tags):
    #     pass

    
    # def rr_around(rr, start=10, stop=10):
    #     rr_0 = np.array(rr_0)
    #     rr_0 = rr_0.astype('float')
    #     rr_1 = np.array(rr_0)
    #     rr_1 = rr_1.astype('float')
    #     rr_1[0:1] = np.nan
    #     rr_2 = np.array(rr_0)
    #     rr_2 = rr_2.astype('float')
    #     rr_2[0:2] = np.nan
    #     rr_3 = np.array(rr_0)
    #     rr_3 = rr_3.astype('float')
    #     rr_3[0:3] = np.nan
    #     rr_4 = np.array(rr_0)
    #     rr_4 = rr_4.astype('float')
    #     rr_4[0:4] = np.nan
    #     rr_5 = np.array(rr_0)
    #     rr_5 = rr_5.astype('float')
    #     rr_5[0:5] = np.nan
    def import_meta(self, path): 
           '''
        This method imports meta information for classyfing, e.g. EEG correlation from a csv file
        '''

    def write_feature_table(self, path):
        '''
        This method writes the feature_table to a file, e.g. a csv file
        '''

    def load_feature_table(self, path):
        '''
        This method loads a feature_table from a file, e.g. a csv file
        '''
        
    def clean_feature_table(self):
        '''
        This method cleans the feature table from NaNs
        '''
        self.features.dropna()

    def select_features(self, feature_list):
        '''
        This method drops all not needed features from the feature_table
        '''

    def evaluate_features(self, feature_list):
        '''
        This method checks the correlation of selected features
        '''

'''
This is the high_level_features.py file
'''

import pandas as pd
from medipy.interfaces.features import Features

class HighLevelFeatures(Features):
    '''
    This is the high_level_features class
    '''
    def __init__(self):
        self.features = pd.DataFrame()

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

'''
This is the features.py file
'''
import abc


class Features(metaclass=abc.ABCMeta):
    '''
    This is the features interface
    '''

    @abc.abstractmethod
    def import_meta(self, path):
        '''
        This method imports meta information for classyfing, e.g. EEG correlation from a csv file
        '''

        raise NotImplementedError

    @abc.abstractmethod
    def write_feature_table(self, path):
        '''
        This method writes the feature_table to a file, e.g. a csv file
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def load_feature_table(self, path):
        '''
        This method loads a feature_table from a file, e.g. a csv file
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def clean_feature_table(self):
        '''
        This method cleans the feature table from NaNs
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def select_features(self, feature_list):
        '''
        This method drops all not needed features from the feature_table
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def evaluate_features(self, feature_list):
        '''
        This method checks the correlation of selected features
        '''
        raise NotImplementedError



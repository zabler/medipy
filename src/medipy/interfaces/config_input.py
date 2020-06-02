'''
This is the config_input.py file
'''
import abc


class ConfigInput(metaclass=abc.ABCMeta):
    '''
    This is the ConfigInput interface
    '''

    @abc.abstractmethod
    def create_config(self):
        '''
        returns object from type Config containing parsed and validated input information
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def validate(self):
        '''
        returns true if the input parameters are valid
        '''
        raise NotImplementedError

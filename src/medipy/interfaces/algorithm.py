'''
This is the algortihm.py file
'''
import abc


class Algorithm(metaclass=abc.ABCMeta):
    '''
    This is the algortihm interface
    '''

    @abc.abstractmethod
    def add_channel_data(self, data):
        '''
        adds only channels to the algorithm, which are required
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def process_data(self):
        '''
        runs the algorithm for analyzing the data if buffer has correct size of filled data
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def plot(self):
        '''
        plots the data in an extra plot
        '''
        raise NotImplementedError

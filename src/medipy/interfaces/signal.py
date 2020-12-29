'''
signal.py
'''
import abc


class Signal(metaclass=abc.ABCMeta):
    '''
    This is the signal interface
    '''
    @abc.abstractmethod
    def plot_signal(self):
        '''
        plots a part of the signal in an extra plot
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def save_signal(self, file):
        '''
        saves the signal object as pickle
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def load_signal(file):
        '''
        loads the signal object from a pickle
        '''
        raise NotImplementedError
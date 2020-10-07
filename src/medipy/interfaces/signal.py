'''
This is the signal.py file
'''
import abc


class Signal(metaclass=abc.ABCMeta):
    '''
    This is the signal interface
    '''
    # Abstract Initation of Attributes?? 

    @abc.abstractmethod
    def plot_signal(self):
        '''
        plots a part of the signal in an extra plot
        '''
        raise NotImplementedError

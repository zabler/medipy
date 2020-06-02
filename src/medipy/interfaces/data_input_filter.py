'''
This is the data_input_filter.py file
'''
import abc


class DataInputFilter(metaclass=abc.ABCMeta):
    '''
    This is the DataInputFilter interface
    '''

    @abc.abstractmethod
    def on_data_read(self, data):
        '''
        runs when DataInput passes data
        preprocesses and filters the data
        returns void
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def on_data_input_finished(self):
        '''
        runs when DataInput has finished
        returns void
        '''
        raise NotImplementedError

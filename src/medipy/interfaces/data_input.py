'''
This is the data_input.py file
'''
import abc


class DataInput(metaclass=abc.ABCMeta):
    '''
    This is the DataInput interface
    '''

    @abc.abstractmethod
    def register_data_input_filter(self, data_input_filter):
        '''
        registers objects from type DataInputFilter
        returns void
        '''
        raise NotImplementedError

    @abc.abstractmethod
    def run(self, input_path, config_path):
        '''
        starts streaming process
        reads incoming data
        passes data to DataInputFilters
        returns void
        '''
        raise NotImplementedError

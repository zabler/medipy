'''
This is the data_output.py file
'''
import abc


class DataOutput(metaclass=abc.ABCMeta):
    '''
    This is the DataOutput interface
    '''

    @abc.abstractmethod
    def on_data_processed(self, task_id, processed_data):
        '''
        starts a process that shows the processed data
        returns void
        '''
        raise NotImplementedError

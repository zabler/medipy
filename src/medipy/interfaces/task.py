'''
This is the task.py file
'''
import abc


class Task(metaclass=abc.ABCMeta):
    '''
    This is the Task interface
    '''

    @abc.abstractmethod
    def add_algorithm(self, name, parameters):
        '''
        '''

        raise NotImplementedError

    @abc.abstractmethod
    def register_data_output(self, data_output):
        '''
        '''
        raise NotImplementedError

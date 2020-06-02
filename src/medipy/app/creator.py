'''
This is the creator.py file
'''
from pyda.app.cmd import CMD
from pyda.app.processor import Processor
from pyda.app.stream_input import StreamInput
from pyda.app.task_filter import TaskFilter


class Creator():
    '''
    This is the creator class
    '''

    def create_data_input(self):
        '''
        returns object of type DataInput
        '''
        return StreamInput()

    def create_data_input_filter(self, filter_string):
        '''
        returns object of type DataInputFilter
        '''
        return TaskFilter(filter_string)

    def create_task(self, task_name):
        '''
        returns object of type Task
        '''
        return Processor(task_name)

    def create_data_output(self):
        '''
        returns object of type DataOutput
        '''
        return CMD()

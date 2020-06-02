'''
This is the task_filter.py file
'''
import ast
from pyda.interfaces.data_input_filter import DataInputFilter


class TaskFilter(DataInputFilter):
    '''
    This is the TaskFilter class
    '''

    def __init__(self, filter_string):
        self.filter_string = filter_string
        self.task = None

    def register_task(self, task):
        '''
        registers task from type Task
        returns void
        '''
        self.task = task

    def on_data_read(self, data):
        '''
        verifies if the data is assigned for the filter's task
        croppes the data into the right format
        passes the data to its task by calling task's read function
        returns void
        '''
        if self.filter_string == data[:len(self.filter_string)]:
            crop = data.split('[')[1]
            crop = crop.split(']')[0]
            try:
                crop = list(ast.literal_eval(crop))
            except ValueError:
                crop = crop.split(', ')

            self.task.on_data_read(crop)

    def on_data_input_finished(self):
        '''
        calls task's finish function
        returns void
        '''
        self.task.on_data_input_finished()

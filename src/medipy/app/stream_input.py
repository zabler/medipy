'''
This is the stream_input.py file
'''
import subprocess
from pyda.interfaces.data_input import DataInput


class StreamInput(DataInput):
    '''
    This is the StreamInput class
    '''

    def __init__(self):
        self.data_input_filters = []

    def register_data_input_filter(self, data_input_filter):
        '''
        registers task_filter from type DataInputFilter
        returns void
        '''
        self.data_input_filters.append(data_input_filter)

    def run(self, input_path, config_path):
        '''
        starts bip2X
        calls filter's read function for every new data line
        calls filter's finish function after last line has been read
        returns void
        '''
        stream = subprocess.Popen(
            ["python", "-m", "bip2X", "-i", input_path, "-c", config_path],
            stdout=subprocess.PIPE, universal_newlines=True)
        for line in stream.stdout:
            for data_input_filter in self.data_input_filters:
                data_input_filter.on_data_read(line)
        for data_input_filter in self.data_input_filters:
            data_input_filter.on_data_input_finished()

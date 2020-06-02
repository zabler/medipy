'''
This is the processor.py file
'''
from pyda.interfaces.task import Task
from pyda.interfaces.data_input_filter import DataInputFilter
from pyda.algorithms.single_mean_time import SingleMeanTime
from pyda.algorithms.multiple_moving_mean_time import MultipleMovingMeanTime


class Processor(Task, DataInputFilter):
    '''
    This is the Processor class
    '''

    def __init__(self, task_name):
        self.task_name = task_name
        self.algorithms = []
        self.data_output = None

    def add_algorithm(self, name, parameters):
        '''
        adds the correct algortihm object including its parameters to processor
        returns void
        '''
        algorithms = {"single_mean_time": SingleMeanTime(name, parameters),
                      "multiple_moving_mean_time": MultipleMovingMeanTime(name, parameters)}

        self.algorithms.append(algorithms[name])

    def register_data_output(self, data_output):
        '''
        buffers the data for possible output
        returns void
        '''
        self.data_output = data_output

    def on_data_read(self, data):
        '''
        passes new data to the algorithms,
        makes the algorithms processing the data if their buffer is full
        passes the processed data to the output
        returns void
        '''
        for algorithm in self.algorithms:
            algorithm.add_channel_data(data)
            if algorithm.buffer_size == len(algorithm.buffer):
                self.data_output.on_data_processed(
                    self.task_name + '_' + algorithm.name + '_' + algorithm.ident,
                    algorithm.process_data())

    def on_data_input_finished(self):
        '''
        trys to make the algortihms process the data for the last time, if there is enough data
        runs the plots function of each algorithm
        returns void
        '''
        for algorithm in self.algorithms:
            try:
                self.data_output.on_data_processed(
                    self.task_name + '_' + algorithm.name + '_' + algorithm.ident,
                    algorithm.process_data())
            except TypeError:
                print(f'Not enough values for process {self.task_name}_{algorithm.name}_{algorithm.ident}')
                print('Last buffer could not be processed')
            finally:
                print(f'Process {self.task_name}_{algorithm.name}_{algorithm.ident} finished')
                # Comment the following for preventing the calculation of algortihm plots
                algorithm.plot()

'''
This is the single_mean_time.py file
'''
import numpy as np
import matplotlib.pyplot as plt
from pyda.interfaces.algorithm import Algorithm


class SingleMeanTime(Algorithm):
    '''
    This is the SingleMeanTime class
    '''

    def __init__(self, name, parameter):
        self.name = name
        self.parameter = parameter

        # From parameter
        self.ident = parameter[0]
        self.channels = [int(x) for x in parameter[1]]
        self.window_seconds = int(parameter[2])
        self.sample_rate_hz = int(parameter[3])

        self.buffer_size = self.window_seconds*self.sample_rate_hz
        self.buffer = []
        self.processed = []

    def add_channel_data(self, data):
        '''
        adds only the channels to the algorithm, which are required
        '''
        channel_data = []
        for channel in self.channels:
            channel_data.append(data[channel])
        self.buffer.append(channel_data)

    def process_data(self):
        '''
        runs the algorithm for analyzing the data, if buffer has the correct size of filled data
        returns the result of the processed data
        a TypeError is raised at the end of the stream input, if the buffer may not be full
        '''
        if self.buffer_size == len(self.buffer):
            calc = []
            for item in self.buffer:
                if item[0] != 'nan':
                    calc.append(item[0])
            result = np.mean(calc)
            self.buffer.clear()
            self.processed.append(result)
            return result
        else:
            raise TypeError

    def plot(self):
        '''
        plots the data for a single channel
        returns void
        '''
        plt.figure()
        plt.plot(self.processed)
        plt.title(f'Mean of {self.ident}')
        plt.xlim(left=0)

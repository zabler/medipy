'''
This is the multiple_moving_mean_time.py file
'''
import numpy as np
import matplotlib.pyplot as plt
from pyda.interfaces.algorithm import Algorithm


class MultipleMovingMeanTime(Algorithm):
    '''
    This is the MultipleMovingMeanTime class
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
        adds only channels to the algorithm, which are required
        returns void
        '''
        channel_data = []
        for channel in self.channels:
            channel_data.append(data[channel])
        self.buffer.append(channel_data)

    def process_data(self):
        '''
        runs the algorithm for analyzing the data, if buffer has the correct size of filled data
        returns the result of the processed data.
        at the end of the stream input, the buffer may not be full, then a TypeError is raised.
        '''
        if self.buffer_size == len(self.buffer):
            result = []
            for ind in range(len(self.channels)):
                channel = [item[ind] for item in self.buffer]
                calc = []
                for val in channel:
                    if val != 'nan':
                        calc.append(val)
                result.append(np.mean(calc))
            del self.buffer[-1]
            self.processed.append(result)
            return result
        else:
            raise TypeError

    def plot(self):
        '''
        plots the processed data for multiple channels
        returns void
        '''
        fig, axs = plt.subplots(len(self.channels), ncols=1)
        fig.suptitle(f'Moving mean of {self.ident}')
        for ind in range(len(self.channels)):
            res = [i[ind] for i in self.processed]
            axs[ind].plot(res)
            axs[ind].set_xlim(left=0)

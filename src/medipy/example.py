'''
This is the __main__.py file
'''
import time
import matplotlib.pyplot as plt
from pyda.app.argument_parser import ArgumentParser
from pyda.app.creator import Creator


class Application():
    '''
    This is the Application class
    '''

    def __init__(self, config):
        self.config = config
        self.stream_input = None
        self.task_filters = []
        self.processors = []
        self.output = None

    def set_up_components(self):
        '''
        creates all necessary components
        1 data_input for all
        k data_input_filters (1 for each stream)
        k tasks (1 for each stream)
        m algortihms (1 for each algorithm of a task)
        1 data_outpout for all
        returns void
        '''
        self.stream_input = Creator().create_data_input()

        for stream in self.config.json_stream_tasks:
            self.task_filters.append(Creator().create_data_input_filter(stream["name"]))
            if stream["algorithms"] is not None:
                self.processors.append(Creator().create_task(stream["name"]))
                for algorithm in stream["algorithms"]:
                    self.processors[-1].add_algorithm(algorithm["name"], algorithm["parameters"])

        self.output = Creator().create_data_output()

    def registrate_components(self):
        '''
        introduces the different components to each other
        returns void
        '''
        # Registration between stream_input and task_filters
        for task_filter in self.task_filters:
            self.stream_input.register_data_input_filter(task_filter)

        # Registration between task filters and corresponding processors
        for ind, task_filter in enumerate(self.task_filters):
            task_filter.register_task(self.processors[ind])

        # Registration between processors and data_output
        for processor in self.processors:
            processor.register_data_output(self.output)

    def run(self):
        '''
        main loop of the application
        returns void
        '''
        self.stream_input.run(self.config.input_path, self.config.config_path)


if __name__ == "__main__":
    START_TIME = time.time()

    CONFIG = ArgumentParser().create_config()
    APP = Application(CONFIG)
    APP.set_up_components()
    APP.registrate_components()
    APP.run()

    ELAPSED_TIME = time.time() - START_TIME
    print("\nConversion finished in " + str(ELAPSED_TIME) + " seconds")

    # Comment following line for hidding single plots of algorithms
    plt.show()

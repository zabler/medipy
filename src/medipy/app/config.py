'''
This is the config.py file
'''


class Config():
    '''
    This is the Config class
    '''

    def __init__(self, args, json_file):
        self.input_path = args.input
        self.config_path = args.config
        self.json_stream_tasks = []
        for task in json_file['tasks']:
            if (task['export']) == "stream":
                self.json_stream_tasks.append(task)

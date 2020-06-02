'''
This is the cmd.py file
'''
from pyda.interfaces.data_output import DataOutput


class CMD(DataOutput):
    '''
    This is the cmd class
    '''

    def on_data_processed(self, task_id, processed_data):
        '''
        prints all processed data to stdout including their task ID
        returns void
        '''
        print(task_id, processed_data)

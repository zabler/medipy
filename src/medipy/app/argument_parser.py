'''
This is the argument_parser.py file
'''
import argparse
import sys
import json
from pyda.app.config import Config
from pyda.interfaces.config_input import ConfigInput


class ArgumentParser(ConfigInput):
    '''
    This is the ArgumentParser class
    '''

    def __init__(self):
        self.args = None
        self.json_config_file = None

        self.parser = argparse.ArgumentParser(description='')
        self.parser.add_argument('-i', '--input', dest='input', required=False,
                                 default=None, type=str, help='Path to input')
        self.parser.add_argument('-c', '--config', dest='config', required=False,
                                 default=None, type=str, help='Path to JSON confg file')

    def parse_cmd(self):
        '''
        parses the cmd line input
        returns void
        '''
        self.args = self.parser.parse_args()

    def load_json(self):
        '''
        loads JSON config file
        raises an ValueError, if JSON file is not valid
        returns void
        '''
        with open(self.args.config) as json_file:
            try:
                self.json_config_file = json.load(json_file)
            except ValueError:
                sys.exit("Config file contains no valid JSON")

    def validate(self):
        '''
        returns true if the input parameters of cmd line are valid
        '''
        if not hasattr(self.args, 'config'):
            sys.exit('No config path found')

        if not hasattr(self.args, 'input'):
            sys.exit('No input path found')

        return True

    def create_config(self):
        '''
        returns object from type Config containing parsed and validated input information from cmd line and JSON
        '''
        self.parse_cmd()
        self.load_json()
        if self.validate():
            return Config(self.args, self.json_config_file)
        sys.exit('Input is not valid')

import json
import os
import re

# try:
#     from const_vals import *
#
# except ModuleNotFoundError:
#     from sharedlib.const_vals import *

JKEY = 'obj'


def get_json_data(file_name, key=None):
    """
    This function uses the python JSON library in order to parse JSON files into usable python objects. Can return
    lists or dictionaries, depending on the JSON file's structure.
    """
    if key is None: key = JKEY
    data = ''
    try:
        with open(file_name) as jsonFile:
            data = json.load(jsonFile)
    except FileNotFoundError:
        file_name = '..' + os.sep + file_name
        with open(file_name) as jsonFile:
            data = json.load(jsonFile)
    return data[key]

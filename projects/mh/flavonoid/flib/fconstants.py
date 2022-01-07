# custom project library imports
from flib.futil import *
from flib.data_types import Flav

# other imports
import os
import sys

SEP = os.sep  # get the right slash. / for linux & mac, \ for windows

# misc strings
CSV = '.csv'
GDATA = 'Gene_data_'
JKEY = 'obj'
NIX = ''
NL = '\n'
SP = ' '
TF = '.txt'
URL_DBGET = 'https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+'

# output directories
DIR_CHEM = SEP + 'Chemical_Data'
DIR_FASTA = SEP + 'FASTA_Data'
DIR_GENE = SEP + 'Gene_Data'
DIR_RGENE = DIR_GENE + SEP + 'Raw_Gene'
DIR_RFASTA = DIR_FASTA + SEP + 'Raw_FASTA'

# pathways for the json data  (and the name of the json object used in all files)
JSON_DIR = '..' + SEP + 'json_data' + SEP
DIR_FJSON = 'fjson' + SEP

JSON_FLAVS = DIR_FJSON + 'flav_names.json'
JSON_FLAV_REL = DIR_FJSON + 'flav_related.json'
JSON_FLAV_SYN = DIR_FJSON + 'flav_syns.json'
JSON_PATH_DICT = DIR_FJSON + 'path_codes_names.json'
JSON_PATH_LIST = DIR_FJSON + 'path_codes.json'
JSON_PLANT_DICT = DIR_FJSON + 'plant_names_codes.json'
JSON_PLANT_DICT_COMMON = DIR_FJSON + 'plant_names_codes_plus.json'
JSON_PLANT_LIST = DIR_FJSON + 'plant_codes.json'
JSON_TST_MED = DIR_FJSON + 'test_med.json'
JSON_TST_SHORT = DIR_FJSON + 'test_short.json'
JSON_TST_SINGLE = DIR_FJSON + 'test_single.json'
JSON_FLAV_NAMES = DIR_FJSON + 'flav_info.json'

# lists and dictionaries made from JSON files
flav_list = get_json_data(JSON_FLAVS, JKEY)  # list of flavonoids of interest (FOI)
flav_names_info = get_json_data(JSON_FLAV_NAMES)
flav_relatives = get_json_data(JSON_FLAV_REL)  # relatives of FOI
flav_synonyms = get_json_data(JSON_FLAV_SYN)  # synonyms for FOI
path_map_dict = get_json_data(JSON_PATH_DICT, JKEY)  # pathway names and codes
path_map_list = get_json_data(JSON_PATH_LIST, JKEY)  # pathway codes
plant_dict_common = get_json_data(JSON_PLANT_DICT_COMMON, JKEY)  # common & sci. names
plant_dict_reg = get_json_data(JSON_PLANT_DICT, JKEY)  # scientific names and plant codes
plant_full_list = get_json_data(JSON_PLANT_LIST, JKEY)  # full list of plant codes
test_med = get_json_data(JSON_TST_MED, JKEY)  # testing list of plant codes
test_short = get_json_data(JSON_TST_SHORT, JKEY)  # testing list of plant codes
test_single = get_json_data(JSON_TST_SINGLE, JKEY)  # testing list of plant codes

plant_dict = plant_dict_reg  # variable exists for ease of value change when testing
plant_list = plant_full_list  # variable exists for ease of value change when testing

flav_data_lists = []
for k in flav_names_info:
    info = flav_names_info.get(k)
    tmp_flav = Flav(k, [], info.get('file'), info.get('code'))
    flav_data_lists.append(tmp_flav)

# keys for accessing dictionaries
EKY = 'EC'
GKY = 'GENE'
NKY = 'NTSEQ'
OKY = 'ORTHOLOGY'
PKY = 'PLANT'

README = SEP + 'ReadMe.txt'

# EC numbers
E01 = 'EC:4.3.1.24'
E02 = 'EC:4.3.1.25'
E03 = 'EC:1.14.14.91'
E04 = 'EC:6.2.1.12'
E05 = 'EC:2.3.1.170'
E06 = 'EC:2.3.1.133'
E07 = 'EC:1.14.14.96'
E08 = 'EC:1.14.13.-'
E09 = 'EC:2.3.1.74'
E10 = 'EC:5.5.1.6'
E11 = 'EC:1.14.20.5'
E12 = 'EC:1.14.19.76'
E13 = 'EC:1.14.14.81'
E14 = 'EC:1.14.14.82'
E15 = 'EC:1.14.11.9'
E16 = 'EC:1.14.20.6'
E17_1 = 'EC:1.1.1.219'
E17_2 = 'EC:1.1.1.234'
E17_FULL = 'EC:1.1.1.219 1.1.1.234'
E18 = 'EC:1.14.20.4'
E19 = 'EC:1.3.1.77'
E20 = 'EC:1.17.1.3'
E21 = 'EC:1.14.14.87'
E22 = 'EC:4.2.1.105'
E23 = 'EC:2.4.1.357'
E24 = 'EC:1.3.1.117'
E_DEC = 'EC:2.3.1.70'
E_GGT = 'EC:2.4.1.74'
E_SOA = 'EC:2.3.1.30'
E_V1G = 'EC:2.4.1.136'

def gfo(lbl, lst):
    for item in lst:
        if lbl == item.code: return item

import os
import sys
from flib.futil import *

SEP = os.sep  # get the right slash. / for linux, \ for windows
# misc strings
NL, SP, NIX, CSV, GDATA = '\n', ' ', '', '.csv', 'Gene_data_'
# output directories
CHEM_DIR, FASTA_DIR, GENE_DIR = SEP + 'Chemical_Data', SEP + 'FASTA_Data', SEP + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
JKEY = 'obj'
JSON_DIR = '..' + SEP + 'json_data' + SEP
F_JSON_DIR = '..' + SEP + 'flavonoid' + SEP + 'fjson' + SEP

FN_FLAVS = F_JSON_DIR + 'flav_names.json'
FN_FLAV_REL = F_JSON_DIR + 'flav_related.json'
FN_FLAV_SYN = F_JSON_DIR + 'flav_syns.json'
FN_PATH_DICT = F_JSON_DIR + 'F_path_names_and_codes.json'
FN_PATH_LIST = F_JSON_DIR + 'F_path_codes.json'
FN_PLANT_DICT = F_JSON_DIR + 'plant_names_codes.json'
FN_PLANT_DICT_COMMON = F_JSON_DIR + 'F_org_names_and_codes_xtra.json'
FN_PLANT_LIST = F_JSON_DIR + 'F_org_codes.json'
FN_TST_MED = F_JSON_DIR + 'F_test_med.json'
FN_TST_SHORT = F_JSON_DIR + 'F_test_short.json'
FN_TST_SINGLE = F_JSON_DIR + 'F_test_single.json'
# keys for accessing dictionaries
O_KEY, E_KEY, N_KEY, P_KEY, G_KEY = 'ORTHOLOGY', 'EC', 'NTSEQ', 'PLANT', 'GENE'
# output file names
FN_AGI = 'apigenin.csv'
FN_BUN = 'butein.csv'
FN_DEC = 'EC-2-3-1-70.csv'
FN_EC = 'epicatechin.csv'
FN_EGT = 'epigallocatechin.csv'
FN_ERD = 'eriodictyol.csv'
FN_GC = 'gallocatechin.csv'
FN_GEN = 'genistein.csv'
FN_GGT = 'EC-2-4-1-74.csv'
FN_HCC = 'isoliquiritigenin.csv'
FN_HWB = 'cyanidin.csv'
FN_KMP = 'kaempferol.csv'
FN_KXN = 'catechin.csv'
FN_LU2 = 'luteolin.csv'
FN_MYC = 'myricetin.csv'
FN_NAR = 'naringenin.csv'
FN_QUER = 'quercetin.csv'
FN_README = SEP + 'ReadMe.txt'
FN_SOA = 'EC-2-3-1-30.csv'
FN_V1G = 'EC-2-4-1-136.csv'

# labels for each of the compounds
AGI = 'Apigenin'
BUN = 'Butein'
CAQ = 'Catechol'
EC = 'Epicatechin'
EGT = 'Epigallocatechin'
ERD = 'Eriodictyol'
GC = 'Gallocatechin'
GEN = 'Genistein'
HCC = 'Isoliquiritigenin'
HWB = 'Cyanidin'
KMP = 'Kaempferol'
KXN = 'Catechin'
LU2 = 'Luteolin'
MYC = 'Myricetin'
NAR = 'Naringenin'
PYG = 'Pyrogallol'
QUE = 'Quercetin'
RCO = 'Resorcinol'
STL = 'Resveratrol'

E01, E02, E03, E04, E05 = 'EC:4.3.1.24', 'EC:4.3.1.25', 'EC:1.14.14.91', 'EC:6.2.1.12', 'EC:2.3.1.170'
E06, E07, E08, E09, E10 = 'EC:2.3.1.133', 'EC:1.14.14.96', 'EC:1.14.13.-', 'EC:2.3.1.74', 'EC:5.5.1.6'
E11, E12, E13, E14, E15 = 'EC:1.14.20.5', 'EC:1.14.19.76', 'EC:1.14.14.81', 'EC:1.14.14.82', 'EC:1.14.11.9'
E16, E17, E17_2, E18, E19 = 'EC:1.14.20.6', 'EC:1.1.1.219', 'EC:1.1.1.219 1.1.1.234', 'EC:1.14.20.4', 'EC:1.3.1.77'
E20, E21, E22, GGT, DEC = 'EC:1.17.1.3', 'EC:1.14.14.87', 'EC:4.2.1.105', 'EC:2.4.1.74', 'EC:2.3.1.70'
SOA, V1G = 'EC:2.3.1.30', 'EC:2.4.1.136'

flav_list = get_json_data(FN_FLAVS, JKEY)
flav_synonyms = get_json_data(FN_FLAV_SYN)
flav_relatives = get_json_data(FN_FLAV_REL)

# fill lists and dictionaries from json files
test_single = get_json_data(FN_TST_SINGLE, JKEY)
test_short = get_json_data(FN_TST_SHORT, JKEY)
test_med = get_json_data(FN_TST_MED, JKEY)

# species codes for plants of interest
plant_full_list = get_json_data(FN_PLANT_LIST, JKEY)
plant_list = plant_full_list

# Dictionary/list of key-value pairs that defines the appropriate genus species for each code
plant_dict_reg = get_json_data(FN_PLANT_DICT, JKEY)
plant_dict_common = get_json_data(FN_PLANT_DICT_COMMON, JKEY)
plant_dict = plant_dict_reg

# list of pathways of interest
path_map_list = get_json_data(FN_PATH_LIST, JKEY)

# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
path_map_dict = get_json_data(FN_PATH_DICT, JKEY)

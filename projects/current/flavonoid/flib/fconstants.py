import os
import sys
from flib.futil import *

SEP = os.sep  # get the right slash. / for linux & mac, \ for windows

# misc strings
CSV = '.csv'
GDATA = 'Gene_data_'
JKEY = 'obj'
NIX = ''
NL = '\n'
SP = ' '
URL_DBGET = 'https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+'

# output directories
DIR_CHEM = SEP + 'Chemical_Data'
DIR_FASTA = SEP + 'FASTA_Data'
DIR_GENE = SEP + 'Gene_Data'
DIR_RAW = DIR_GENE + SEP + 'Raw_Gene'

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

# lists and dictionaries made from JSON files
flav_list = get_json_data(JSON_FLAVS, JKEY)  # list of flavonoids of interest (FOI)
flav_relatives = get_json_data(JSON_FLAV_REL)  # relatives of FOI
flav_synonyms = get_json_data(JSON_FLAV_SYN)  # synonyms for FOI
path_map_dict = get_json_data(JSON_PATH_DICT, JKEY)  # pathway names and codes
path_map_list = get_json_data(JSON_PATH_LIST, JKEY)  # pathway codes
plant_dict_common = get_json_data(JSON_PLANT_DICT_COMMON, JKEY)  # contains plants common and scientific names
plant_dict_reg = get_json_data(JSON_PLANT_DICT, JKEY)  # scientific names and plant codes
plant_full_list = get_json_data(JSON_PLANT_LIST, JKEY)  # full list of plant codes
test_med = get_json_data(JSON_TST_MED, JKEY)  # testing list of plant codes
test_short = get_json_data(JSON_TST_SHORT, JKEY)  # testing list of plant codes
test_single = get_json_data(JSON_TST_SINGLE, JKEY)  # testing list of plant codes

plant_dict = plant_dict_reg  # variable exists for ease of value change when testing
plant_list = plant_full_list  # variable exists for ease of value change when testing

# keys for accessing dictionaries
EKY = 'EC'
GKY = 'GENE'
NKY = 'NTSEQ'
OKY = 'ORTHOLOGY'
PKY = 'PLANT'

# output file names
README = SEP + 'ReadMe.txt'
CSV_AGI = 'apigenin.csv'
CSV_BUN = 'butein.csv'
CSV_DEC = 'EC-2-3-1-70.csv'
CSV_EC = 'epicatechin.csv'
CSV_EGT = 'epigallocatechin.csv'
CSV_ERD = 'eriodictyol.csv'
CSV_GC = 'gallocatechin.csv'
CSV_GEN = 'genistein.csv'
CSV_GGT = 'EC-2-4-1-74.csv'
CSV_HCC = 'isoliquiritigenin.csv'
CSV_HWB = 'cyanidin.csv'
CSV_KMP = 'kaempferol.csv'
CSV_KXN = 'catechin.csv'
CSV_LU2 = 'luteolin.csv'
CSV_MYC = 'myricetin.csv'
CSV_NAR = 'naringenin.csv'
CSV_QUE = 'quercetin.csv'
CSV_SOA = 'EC-2-3-1-30.csv'
CSV_V1G = 'EC-2-4-1-136.csv'

CSV_AZEL = 'afzelechin.csv'
CSV_APIF = 'apiforol.csv'
CSV_BUTN = 'butin.csv'
CSV_CHSN = 'chrysin.csv'
CSV_DDZN = 'daidzein.csv'
CSV_DF74 = "7-4dihydroxyflavone.csv"
CSV_DFV = 'liquiritigenin.csv'
CSV_DHGN = '2Hydroxy_2-3dihydrogenistein.csv'
CSV_DHKM = 'dihydrokaempferol.csv'
CSV_DHMF = 'dihydrotricetin.csv'
CSV_DHMY = 'dihydromyricetin.csv'
CSV_DHQU = 'dihydroquercetin.csv'
CSV_DLM = 'delphinidin.csv'
CSV_EZEL = 'epiafzelechin.csv'
CSV_ERCH = 'eriodictyol_chalcone.csv'
CSV_FSTN = 'dihydrofisetin.csv'
CSV_G50 = 'phloretin.csv'
CSV_GBZL = 'garbanzol.csv'
CSV_GLGN = 'galangin.csv'
CSV_HDZ6 = '6hydroxydaidzein.csv'
CSV_HESP = 'hesperetin.csv'
CSV_LDLM = 'leucodelphinidin.csv'
CSV_LHWB = 'leucocyanidin.csv'
CSV_LPLR = 'leucopelargonidin.csv'
CSV_LUTF = 'luteoforol.csv'
CSV_MYF = 'tricetin.csv'
CSV_NACH = 'naringenin_chalcone.csv'
CSV_PBAN = 'pinobanksin.csv'
CSV_PCCH = 'pinocembrin_chalcone.csv'
CSV_PCEM = 'pinocembrin.csv'
CSV_PLRG = 'pelargonidin.csv'
CSV_T2674 = "2-6-7-4tetrahydroxyisoflavanone.csv"
CSV_T274 = "2-7-4trihydroxyisoflavanone.csv"
CSV_T674 = "6-7-4Trihydroxyflavanone.csv"

# labels for each of the compounds


AGI = 'Apigenin'
APIF = 'Apiforol'
BUN = 'Butein'
BUTN = 'Butin'
CAQ = 'Catechol'
CHSN = 'Chrysin'
DDZN = 'Daidzein'
DF74 = "7,4'-Dihydroxyflavone"
DFV = 'Liquiritigenin'
DHGN = '2-Hydroxy-2,3-dihydrogenistein'
DHKM = 'Dihydrokaempferol'
DHMY = 'Dihydromyricetin'
DHQU = 'Dihydroquercetin'
DLM = 'Delphinidin'
EC = 'Epicatechin'
EGT = 'Epigallocatechin'
ERCT = 'Eriocitrin'
ERD = 'Eriodictyol'
GC = 'Gallocatechin'
GEN = 'Genistein'
HCC = 'Isoliquiritigenin'
HWB = 'Cyanidin'
KMP = 'Kaempferol'
KXN = 'Catechin'
LDLM = 'Leucodelphinidin'
LHWB = 'Leucocyanidin'
LU2 = 'Luteolin'
MYC = 'Myricetin'
NACH = 'Naringenin Chalcone'
NAR = 'Naringenin'
PYG = 'Pyrogallol'
QUE = 'Quercetin'
RCO = 'Resorcinol'
STL = 'Resveratrol'
AZEL = '(+)-Afzelechin'
DHMF = 'Dihydrotricetin'
EZEL = '(-)-Epiafzelechin'
ERCH = 'Eriodictyol Chalcone'
FSTN = 'Dihydrofisetin'
G50 = 'Phloretin'
GBZL = 'Garbanzol'
GLGN = 'Galangin'
HDZ6 = '6-Hydroxydaidzein'
HESP = 'Hesperetin'
LPLR = 'Leucopelargonidin'
LUTF = 'Luteoforol'
MYF = 'Tricetin'
PBAN = 'Pinobanksin'
PCCH = 'Pinocembrin chalcone'
PCEM = 'Pinocembrin'
PLRG = 'Pelargonidin'
T2674 = "2,6,7,4'-Tetrahydroxyisoflavanone"
T274 = "2,7,4'-Trihydroxyisoflavanone"
T674 = "6,7,4'-Trihydroxyflavanone"

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
E17 = 'EC:1.1.1.219'
E17_2 = 'EC:1.1.1.219 1.1.1.234'
E18 = 'EC:1.14.20.4'
E19 = 'EC:1.3.1.77'
E20 = 'EC:1.17.1.3'
E21 = 'EC:1.14.14.87'
E22 = 'EC:4.2.1.105'
E23 = '2.4.1.357'
E24 = '1.3.1.117'
E25 = '1.1.1.234'

GGT = 'EC:2.4.1.74'
DEC = 'EC:2.3.1.70'
SOA = 'EC:2.3.1.30'
V1G = 'EC:2.4.1.136'

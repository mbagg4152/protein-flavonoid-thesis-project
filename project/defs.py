# common/important strings. holding in vars helps with human error
# chem_path, fasta_path & gene_path are filled after file checks
# using ' instead of " since it is common practice unless needing to use actual ' character
from util import *
import datetime
import os

slash = os.sep  # get the right slash. / for linux, \ for windows

# output file names
apig_file = 'apigenin.txt'
bute_file = 'butein.txt'
cate_file = 'catechin.txt'
cyan_file = 'cyanidin.txt'
ecat_file = 'epicatechin.txt'
epig_file = 'epigallocatechin.txt'
erio_file = 'eriodictyol.txt'
gall_file = 'gallocatechin.txt'
kaem_file = 'kaempferol.txt'
lute_file = 'luteolin.txt'
nari_file = 'naringenin.txt'
quer_file = 'quercetin.txt'
readme = slash + 'ReadMe.txt'

# output directories
chem_dir = slash + 'Chemical_Data'
fasta_dir = slash + 'FASTA_Data'
gene_dir = slash + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
cate_test_json = 'info' + slash + 'testCate.json'
ecat_test_json = 'info' + slash + 'testEcat.json'
erio_test_json = 'info' + slash + 'testErio.json'
f_species_dict_json = 'info' + slash + 'fullSpeciesDict.json'
f_species_list_json = 'info' + slash + 'fullSpeciesList.json'
j_key = 'obj'
nari_test_json = 'info' + slash + 'testNari.json'
path_dict_json = 'info' + slash + 'pathCodeDict.json'
path_list_json = 'info' + slash + 'pathCodes.json'

# fill lists and dictionaries from json files
cate_test = get_json_data(cate_test_json, j_key)
ecat_test = get_json_data(ecat_test_json, j_key)
erio_test = get_json_data(erio_test_json, j_key)
nari_test = get_json_data(nari_test_json, j_key)

# main lists for the program
# species codes for plants of interest,  "mtr" cut due to code errors
full_list = get_json_data(f_species_list_json, j_key)

# Dictionary/list of key-value pairs that defines the appropriate genus species for each code
species_pairs = get_json_data(f_species_dict_json, j_key)

# list of pathways of interest
pathway_list = get_json_data(path_list_json, j_key)

# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
pathway_pairs = get_json_data(path_dict_json, j_key)

# enzyme catalyst variables used in the analysis. these are constants
# prefix C for catalyst. each letter used is not important, just labels
CA = 'EC:4.3.1.24'
CB = 'EC:4.3.1.25'
CC = 'EC:1.14.14.91'
CD = 'EC:6.2.1.12'
CE = 'EC:23.1.170'
CF = 'EC:2.3.1.133'
CG = 'EC:1.14.14.96'
CH = 'EC:1.14.13.-'
CI = 'EC:2.3.1.74'
CJ = 'EC:5.5.1.6'
CK = 'EC:1.14.20.5'
CL = 'EC:1.14.19.76'
CM = 'EC:1.14.14.81'
CN = 'EC:1.14.14.82'
CO = 'EC:1.14.11.9'
CP = 'EC:1.14.20.6'
CQ = 'EC:1.1.12.19'
CR = 'EC:1.14.20.4'
CS = 'EC:1.3.1.77'
CT = 'EC:1.17.1.3'

# logical combination of pathways for each chemical
start = (CA or CB) and CC and CD
butein = start and CE
naringenin = start and CI and CJ
eriodictyol = (naringenin and (CM or CN)) or (CH or (CF and CG) and CI)
apigenin = (CK or CL) and naringenin
luteolin = apigenin and (CM or CN)
kaempferol = naringenin and CO and CP
quercetin = eriodictyol and ((CO and CP) or (kaempferol and (CM or CN)))
catechin = eriodictyol and CO and CQ and CT
gallocatechin = eriodictyol and CO and CM and CQ and CT
epigallocatechin = eriodictyol and CO and CM and CQ and CR and CS
cyanidin = eriodictyol and CO and CQ and CR
epicatechin = cyanidin and CS

# misc. string values
CSV = '.csv'
GDATA = 'Gene_data_'
NIX = ""
NL = '\n'
NOW = datetime.datetime.now()
SP = ' '

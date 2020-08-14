# common/important strings. holding in vars helps with human error
# chem_path, fasta_path & gene_path are filled after file checks
# using ' instead of " since it is common practice unless needing to use actual ' character
import json
import datetime
import os

slash = os.sep  # get the right slash. / for linux, \ for windows
cate_file = 'cate-species.txt'
chem_dir = slash + 'Chemical_Data'
csv_ext = '.csv'
ecat_file = 'ecat-species.txt'
erio_file = 'erio-species.txt'
fasta_dir = slash + 'FASTA_Data'
gene_dir = slash + 'Gene_Data'
gene_prefix = 'Gene_data_'
lute_file = 'lute-species.txt'
nari_file = 'nari-species.txt'
nix = ''
nl = '\n'
readme = slash + 'ReadMe.txt'
sp = ' '
f_species_list_json = 'info' + slash + 'fullSpeciesList.json'
f_species_dict_json = 'info' + slash + 'fullSpeciesDict.json'
path_list_json = 'info' + slash + 'pathCodes.json'
path_dict_json = 'info' + slash + 'pathCodeDict.json'
ecat_test_json = 'info' + slash + 'testEcat.json'
cate_test_json = 'info' + slash + 'testCate.json'
erio_test_json = 'info' + slash + 'testErio.json'
nari_test_json = 'info' + slash + 'testNari.json'
j_key = 'obj'

# logical operations
phen_cinn = ['or', 'EC:4.3.1.24', 'EC:4.3.1.25']
cinn_pcoa = ['and', 'EC:6.2.1.12', 'EC:1.14.14.91']
pcoa_ccoa1 = ['and', 'EC:1.14.13.-']
pcoa_ccoa2 = ['and', 'EC:2.3.1.133', 'EC:1.14.14.96']
pcoa_nari = ['and', 'EC:2.3.1.74', 'EC:5.5.1.6']
nari_erio = ['or', 'EC:1.14.14.81', 'EC:1.14.14.82']
ccoa_erio = ['and', 'EC:2.3.1.74']
erio_leuc = ['and', 'EC:1.14.11.9', 'EC:1.1.12.19']
leuc_cate = ['and', 'EC:1.17.1.3']
leuc_cyan = ['and', 'EC:1.14.20.4']
cyan_epic = ['and', 'EC:1.3.1.77']
erio_lute = ['or', 'EC:1.14.20.5', 'EC:1.14.19.76']
nari_apig = ['or', 'EC:1.14.20.5', 'EC:1.14.19.76']
apig_lute = ['or', 'EC:1.14.14.81', 'EC:1.14.14.82']


def get_json_data(file_name, key):
    with open(file_name) as jsonFile:
        data = json.load(jsonFile)
        # print(data[key])
        return data[key]


# fill lists and dictionaries from json files
cate_test = get_json_data(cate_test_json, j_key)
ecat_test = get_json_data(ecat_test_json, j_key)
nari_test = get_json_data(nari_test_json, j_key)
erio_test = get_json_data(erio_test_json, j_key)
# main lists for the program
# species codes for plants of interest,  "mtr" cut due to code errors
full_list = get_json_data(f_species_list_json, j_key)
# Dictionary/list of key-value pairs that defines the appropriate genus species for each code
species_pairs = get_json_data(f_species_dict_json, j_key)

# list of pathways of interest
pathway_list = get_json_data(path_list_json, j_key)
# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
# "mtr": "Medicago truncatula" cut due to code error
pathway_pairs = get_json_data(path_dict_json, j_key)

now = datetime.datetime.now()

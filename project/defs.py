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
bute_file = 'bute.txt'
apig_file = 'apig.txt'
kaem_file = 'kaem.txt'
quer_file = 'quer.txt'
cyan_file = 'cyan.txt'
epig_file = 'epig.txt'
gall_file = 'gall.txt'
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

# shortened path variables, used for simplification
# the letters are not meaningful
eca = 'EC:4.3.1.24'
ecb = 'EC:4.3.1.25'
ecc = 'EC:1.14.14.91'
ecd = 'EC:6.2.1.12'
ece = 'EC:23.1.170'
ecf = 'EC:2.3.1.133'
ecg = 'EC:1.14.14.96'
ech = 'EC:1.14.13.-'
eci = 'EC:2.3.1.74'
ecj = 'EC:5.5.1.6'
eck = 'EC:1.14.20.5'
ecl = 'EC:1.14.19.76'
ecm = 'EC:1.14.14.81'
ecn = 'EC:1.14.14.82'
eco = 'EC:1.14.11.9'
ecp = 'EC:1.14.20.6'
ecq = 'EC:1.1.12.19'
ecr = 'EC:1.14.20.4'
ecs = 'EC:1.3.1.77'
ect = 'EC:1.17.1.3'

start = (eca or ecb) and (ecc and ecd)
butein = start and ece
naringenin = start and eci and ecj
eriodictyol = (naringenin and (ecm or ecn)) or (ech or (ecf and ecg) and eci)
apigenin = naringenin and (eck or ecl)
lute = apigenin and (ecm or ecn)
kaem = naringenin and eco and ecp
quer = eriodictyol and ((eco and ecp) or (kaem and (ecm or ecn)))
cate = eriodictyol and eco and ecq and ect
gall = eriodictyol and eco and ecm and ecq and ect
epig = eriodictyol and eco and ecm and ecq and ecr and ecs
cyan = eriodictyol and eco and ecq and ecr
ecat = cyan and ecs


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

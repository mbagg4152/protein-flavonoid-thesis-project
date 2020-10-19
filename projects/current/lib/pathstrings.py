import os

SEP = os.sep  # get the right slash. / for linux, \ for windows

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
FN_SOA = 'EC-2-3-1-30.csv'
FN_V1G = 'EC-2-4-1-136.csv'
FN_README = SEP + 'ReadMe.txt'

# output directories
CHEM_DIR = SEP + 'Chemical_Data'
FASTA_DIR = SEP + 'FASTA_Data'
GENE_DIR = SEP + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
JSON_DIR = '..' + SEP + 'json_data' + SEP
FLAV_JSON = JSON_DIR + 'flavonoid' + SEP
PROT_JSON = JSON_DIR + 'protein' + SEP
JKEY = 'obj'
FN_FLAVS = FLAV_JSON + 'flavs.json'
FN_PATH_DICT = FLAV_JSON + 'path_names_codes.json'
FN_PATH_LIST = FLAV_JSON + 'path_codes.json'
FN_PLANT_DICT = FLAV_JSON + 'plant_names_codes.json'
FN_PLANT_DICT_COMMON = FLAV_JSON + 'plant_names_codes_more.json'
FN_PLANT_LIST = FLAV_JSON + 'plant_codes.json'
FN_TEST_MED = FLAV_JSON + 'test_med.json'
FN_TEST_SHORT = FLAV_JSON + 'test_short.json'
FN_TEST_SINGLE = FLAV_JSON + 'test_single.json'

FN_LIGAND = PROT_JSON + 'ligand.json'
FN_LIGAND_CODES = PROT_JSON + 'ligand_codes.json'
FN_PDB_IDS = PROT_JSON + 'pdb_ids.json'
FN_PROTEIN = PROT_JSON + 'protein.json'
FN_PDB_IDS_SHORT = PROT_JSON + 'pdb_ids_short.json'

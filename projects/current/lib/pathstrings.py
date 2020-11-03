import os

SEP = os.sep  # get the right slash. / for linux, \ for windows

# output file names
FN_AGI, FN_BUN, FN_DEC, FN_EC = 'apigenin.csv', 'butein.csv', 'EC-2-3-1-70.csv', 'epicatechin.csv'
FN_EGT, FN_ERD, FN_GC, FN_GEN = 'epigallocatechin.csv', 'eriodictyol.csv', 'gallocatechin.csv', 'genistein.csv'
FN_GGT, FN_HCC, FN_HWB, FN_KMP = 'EC-2-4-1-74.csv', 'isoliquiritigenin.csv', 'cyanidin.csv', 'kaempferol.csv'
FN_KXN, FN_LU2, FN_MYC, FN_NAR = 'catechin.csv', 'luteolin.csv', 'myricetin.csv', 'naringenin.csv'
FN_QUER, FN_SOA, FN_V1G, FN_README = 'quercetin.csv', 'EC-2-3-1-30.csv', 'EC-2-4-1-136.csv', SEP + 'ReadMe.txt'

# output directories
CHEM_DIR, FASTA_DIR, GENE_DIR = SEP + 'Chemical_Data', SEP + 'FASTA_Data', SEP + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
JSON_DIR = '..' + SEP + 'json_data' + SEP
FLAV_JSON, PROT_JSON, JKEY = JSON_DIR + 'flavonoid' + SEP, JSON_DIR + 'protein' + SEP, 'obj'
FN_FLAVS, FN_PATH_DICT = FLAV_JSON + 'flavs.json', FLAV_JSON + 'path_names_codes.json'
FN_PATH_LIST, FN_PLANT_DICT = FLAV_JSON + 'path_codes.json', FLAV_JSON + 'plant_names_codes.json'
FN_PLANT_DICT_COMMON, FN_PLANT_LIST = FLAV_JSON + 'plant_names_codes_more.json', FLAV_JSON + 'plant_codes.json'
FN_MED, FN_SHORT, FN_SINGLE = FLAV_JSON + 'test_med.json', FLAV_JSON + 'test_short.json', FLAV_JSON + 'test_single.json'

FN_LIGAND, FN_LIGAND_CODES = PROT_JSON + 'ligand.json', PROT_JSON + 'ligand_codes.json'
FN_PDB_IDS, FN_PDB_IDS_SHORT = PROT_JSON + 'pdb_ids.json', PROT_JSON + 'pdb_ids_short.json'
FN_PROTEIN, FN_PDB_SINGLE = PROT_JSON + 'protein.json', PROT_JSON + 'pdb_id_single.json'

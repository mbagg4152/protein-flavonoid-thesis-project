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
JKEY = 'obj'

JSON_DIR = '..' + SEP + 'json_data' + SEP
JDIR_FLAV = JSON_DIR + 'flavonoid' + SEP
JDIR_PROT = JSON_DIR + 'protein' + SEP
JDIR_NPASS = JDIR_FLAV + 'npass' + SEP
JDIR_TST = JDIR_FLAV + 'test' + SEP

FN_FLAVS = JDIR_FLAV + 'flavs.json'
FN_FLAV_SYN = JDIR_FLAV + 'flav_synonyms.json'
FN_FLAV_REL = JDIR_FLAV + 'flav_relatives.json'
FN_PATH_DICT = JDIR_FLAV + 'path_names_codes.json'
FN_PATH_LIST = JDIR_FLAV + 'path_codes.json'
FN_PLANT_DICT = JDIR_FLAV + 'plant_names_codes.json'
FN_PLANT_DICT_COMMON = JDIR_FLAV + 'plant_names_codes_more.json'
FN_PLANT_LIST = JDIR_FLAV + 'plant_codes.json'
FN_MED = JDIR_TST + 'test_med.json'
FN_SHORT = JDIR_TST + 'test_short.json'
FN_SINGLE = JDIR_TST + 'test_single.json'

FN_LIGAND = JDIR_PROT + 'ligand.json'
FN_LIGAND_CODES = JDIR_PROT + 'ligand_codes.json'
FN_PDB_IDS = JDIR_PROT + 'pdb_ids.json'
FN_PDB_IDS_SHORT = JDIR_PROT + 'pdb_ids_short.json'
FN_PROTEIN = JDIR_PROT + 'protein.json'
FN_PDB_SINGLE = JDIR_PROT + 'pdb_id_single.json'
FN_PLANES = JDIR_PROT + 'planes.json'

FN_NPLANT_DICT = JDIR_NPASS + 'nplant_dict.json'
FN_NPLANT_CODES = JDIR_NPASS + 'nplant_codes.json'

FN_NCOMP_DICT = JDIR_NPASS + 'ncomp_dict.json'
FN_NCOMP_CODES = JDIR_NPASS + 'ncomp_codes.json'

FN_SMILES = JDIR_PROT + 'smiles.json'
FN_SMILES_PART = JDIR_PROT + 'smiles_part.json'
FN_KNAP = JDIR_PROT + 'knap_ids.json'

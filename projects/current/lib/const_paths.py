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
FN_README = SEP + 'ReadMe.txt'
FN_SOA = 'EC-2-3-1-30.csv'
FN_V1G = 'EC-2-4-1-136.csv'

# output directories
CHEM_DIR, FASTA_DIR, GENE_DIR = SEP + 'Chemical_Data', SEP + 'FASTA_Data', SEP + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
JKEY = 'obj'
JSON_DIR = '..' + SEP + 'json_data' + SEP

FN_FLAVS = JSON_DIR + 'F_flav_names.json'
FN_FLAV_REL = JSON_DIR + 'F_flav_relatives.json'
FN_FLAV_SYN = JSON_DIR + 'F_flav_synonyms.json'
FN_PATH_DICT = JSON_DIR + 'F_path_names_and_codes.json'
FN_PATH_LIST = JSON_DIR + 'F_path_codes.json'
FN_PLANT_DICT = JSON_DIR + 'F_org_names_and_codes.json'
FN_PLANT_DICT_COMMON = JSON_DIR + 'F_org_names_and_codes_xtra.json'
FN_PLANT_LIST = JSON_DIR + 'F_org_codes.json'
FN_TST_MED = JSON_DIR + 'F_test_med.json'
FN_TST_SHORT = JSON_DIR + 'F_test_short.json'
FN_TST_SINGLE = JSON_DIR + 'F_test_single.json'

FN_LIGAND_CODES = JSON_DIR + 'P_ligand_codes.json'
FN_LIGAND_INFO = JSON_DIR + 'P_ligand_info.json'
FN_LIG_IDENTIFIERS = JSON_DIR + 'P_ligand_identifiers.json'
FN_LIG_IDENTIFIERS_PART = JSON_DIR + 'P_ligand_identifiers_part.json'
FN_PDB_IDS = JSON_DIR + 'P_pdb_ids.json'
FN_PDB_IDS_SHORT = JSON_DIR + 'P_pdb_ids_short.json'
FN_PDB_KNAP_IDS = JSON_DIR + 'P_knapsack_ids.json'
FN_PDB_SINGLE = JSON_DIR + 'P_pdb_single_id.json'
FN_PROTEIN = JSON_DIR + 'P_protein_ids.json'
FN_RINGS = JSON_DIR + 'P_ligand_rings.json'
FN_LIG_TESTS = JSON_DIR + 'P_test_ligand_ids.json'
FN_PUBCHEM_PDB = JSON_DIR + 'P_pubchem_pdb.json'

FN_NCOMP_CODES = JSON_DIR + 'N_comp_codes.json'
FN_NCOMP_DICT = JSON_DIR + 'N_comp_codes_and_names.json'
FN_NPLANT_CODES = JSON_DIR + 'N_org_codes.json'
FN_NPLANT_DICT = JSON_DIR + 'N_org_codes_and_names.json'

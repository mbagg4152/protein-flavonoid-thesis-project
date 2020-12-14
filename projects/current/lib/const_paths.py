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

FN_FLAVS = JSON_DIR + 'FlavNames.json'
FN_FLAV_REL = JSON_DIR + 'FlavRelatives.json'
FN_FLAV_SYN = JSON_DIR + 'FlavSynonyms.json'
FN_PATH_DICT = JSON_DIR + 'FlavPathNamesAndCodes.json'
FN_PATH_LIST = JSON_DIR + 'FlavPathCodes.json'
FN_PLANT_DICT = JSON_DIR + 'FlavOrgNamesAndCodes.json'
FN_PLANT_DICT_COMMON = JSON_DIR + 'FlavOrgNamesAndCodesExtra.json'
FN_PLANT_LIST = JSON_DIR + 'FlavPlantCodes.json'
FN_TST_MED = JSON_DIR + 'FlavTestMed.json'
FN_TST_SHORT = JSON_DIR + 'FlavTestShort.json'
FN_TST_SINGLE = JSON_DIR + 'FlavTestSingle.json'

FN_LIGAND_CODES = JSON_DIR + 'ProtLigandCodes.json'
FN_LIGAND_INFO = JSON_DIR + 'ProtLigandInfo.json'
FN_LIG_IDENTIFIERS = JSON_DIR + 'ProtLigandIdentifiers.json'
FN_LIG_IDENTIFIERS_PART = JSON_DIR + 'ProtLigandIdentifiersPartial.json'
FN_PDB_IDS = JSON_DIR + 'ProtPdbIDs.json'
FN_PDB_IDS_SHORT = JSON_DIR + 'ProtPdbIDsShort.json'
FN_PDB_KNAP_IDS = JSON_DIR + 'ProtKnapsackIDs.json'
FN_PDB_SINGLE = JSON_DIR + 'ProtPdbSingleID.json'
FN_PROTEIN = JSON_DIR + 'ProtProteinIDs.json'
FN_RINGS = JSON_DIR + 'ProtRings.json'
FN_LIG_TESTS = JSON_DIR + 'ProtTestLigandIDs.json'

FN_NCOMP_CODES = JSON_DIR + 'NpassCompCodes.json'
FN_NCOMP_DICT = JSON_DIR + 'NpassCompCodesNames.json'
FN_NPLANT_CODES = JSON_DIR + 'NpassOrgCodes.json'
FN_NPLANT_DICT = JSON_DIR + 'NpassOrgCodesNames.json'

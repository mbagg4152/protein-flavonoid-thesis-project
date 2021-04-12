import os

SEP = os.sep
JKEY = 'obj'
P_JSON_DIR = '..' + SEP + 'protein' + SEP + 'pjson' + SEP
FN_LIGAND_CODES = P_JSON_DIR + 'P_ligand_codes.json'
FN_LIGAND_INFO = P_JSON_DIR + 'P_ligand_info.json'
FN_LIG_IDENTIFIERS = P_JSON_DIR + 'P_ligand_identifiers.json'
FN_LIG_IDENTIFIERS_PART = P_JSON_DIR + 'P_ligand_identifiers_part.json'
FN_PDB_IDS = P_JSON_DIR + 'P_pdb_ids.json'
FN_PDB_IDS_SHORT = P_JSON_DIR + 'P_pdb_ids_short.json'
FN_PDB_KNAP_IDS = P_JSON_DIR + 'P_knapsack_ids.json'
FN_PDB_SINGLE = P_JSON_DIR + 'P_pdb_single_id.json'
FN_PROTEIN = P_JSON_DIR + 'P_protein_ids.json'
FN_RINGS = P_JSON_DIR + 'P_ligand_rings.json'
FN_LIG_TESTS = P_JSON_DIR + 'P_test_ligand_ids.json'
FN_PUBCHEM_PDB = P_JSON_DIR + 'P_pubchem_pdb.json'

CHUNK_SIZE = 6
PART_URL = "https://files.rcsb.org/view/"

# PDB file keys
K_ATM = 'ATOM'
K_CMP = 'COMPND'
K_EC = 'EC:'
K_EX_SYS = 'EXPRESSION_SYSTEM'
K_EX_TAX = 'EXPRESSION_SYSTEM_TAXID'
K_HAT = 'HETATM'
K_HEAD = 'HEADER'
K_ORG = 'ORGAN'
K_ORG_CMN = 'ORGANISM_COMMON'
K_ORG_SCI = 'ORGANISM_SCIENTIFIC'
K_ORG_TAX = 'ORGANISM_TAXID'
K_REV = 'REVDAT'
K_SRC = 'SOURCE'
K_TTL = 'TITLE'

# regular expressions
RE_EC = r'([0-9]\.{1}[^ ,;]*)'  # look for EC number format
RE_XTRA_SP = r' {2,}'  # look for 2+ spaces
RE_REC_VAL = r':(.*);'  # record values are between : and ; for specific records
RE_WORDS = r'\s*(\S[\S| ]*\S)\s*'

# directory/path values

path_cwd = os.getcwd()
out_dir = path_cwd + SEP + 'pdb-output' + SEP
image_dir = path_cwd + SEP + 'images' + SEP
formatted_out = out_dir + 'format_out_all.txt'
formatted_out_small = out_dir + 'format_out_partial.txt'
formatted_basic = out_dir + 'basic_info_formatted.txt'
pdb_dir = out_dir + 'pdb_files' + SEP
pdb_test_dir = path_cwd + SEP + 'pdb-test-files' + SEP

# default values
DS = ' '  # default string
DF = 0.0  # default float
DI = 0  # default int
DL = []  # default list
NS = 'NONE'

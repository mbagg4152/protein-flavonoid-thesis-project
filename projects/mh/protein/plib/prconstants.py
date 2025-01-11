import os
from plib.pdb_util import *

SEP = os.sep
JKEY = 'obj'
P_JSON_DIR = '..' + SEP + 'protein' + SEP + 'prjson' + SEP
FN_LIGAND_CODES = P_JSON_DIR + 'lig_codes.json'
FN_LIGAND_INFO = P_JSON_DIR + 'lig_info.json'
FN_LIG_IDENTIFIERS = P_JSON_DIR + 'lig_identifiers.json'
FN_LIG_IDENTIFIERS_PART = P_JSON_DIR + 'lig_identifiers_part.json'
FN_PDB_IDS = P_JSON_DIR + 'pdb_ids.json'
FN_PDB_IDS_SHORT = P_JSON_DIR + 'pdb_ids_short.json'
FN_PDB_KNAP_IDS = P_JSON_DIR + 'knap_ids.json'
FN_PDB_SINGLE = P_JSON_DIR + 'pdb_single_id.json'
FN_PROTEIN = P_JSON_DIR + 'protein_ids.json'
FN_RINGS = P_JSON_DIR + 'lig_rings.json'
FN_LIG_TESTS = P_JSON_DIR + 'test_lig_ids.json'
FN_PUBCHEM_PDB = P_JSON_DIR + 'pubchem_pdb.json'
FN_TESTNAMES = P_JSON_DIR + 'test_pdb_filenames.json'

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
pdb_test_dir = path_cwd + SEP + 'local-pdb-files' + SEP
pdb_m0_only = pdb_test_dir + 'test_m0_files' + SEP

pdb_filenames = get_json_data(FN_TESTNAMES)
ligand_info = get_json_data(FN_LIGAND_INFO)
ligand_codes = get_json_data(FN_LIGAND_CODES)
pdb_id_list = get_json_data(FN_PDB_IDS, JKEY)
pdb_id_list_short = get_json_data(FN_PDB_IDS_SHORT)
pdb_id_single = get_json_data(FN_PDB_SINGLE)
protein_info = get_json_data(FN_PROTEIN)
struct_rings = get_json_data(FN_RINGS)

# default values
DS = ' '  # default string
DF = 0.0  # default float
DI = 0  # default int
DL = []  # default list
NS = 'NONE'

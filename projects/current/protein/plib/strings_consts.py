import os

CHUNK_SIZE = 6
PART_URL = "https://files.rcsb.org/view/"

# PDB file keys
K_CMP = 'COMPND'
K_EC = 'EC:'
K_EX_SYS = 'EXPRESSION_SYSTEM'
K_EX_TAX = 'EXPRESSION_SYSTEM_TAXID'
K_HEAD = 'HEADER'
K_REV = 'REVDAT'
K_ORG = 'ORGAN'
K_ORG_CMN = 'ORGANISM_COMMON'
K_ORG_SCI = 'ORGANISM_SCIENTIFIC'
K_ORG_TAX = 'ORGANISM_TAXID'
K_SRC = 'SOURCE'
K_TTL = 'TITLE'
K_ATM = 'ATOM'
K_HAT = 'HETATM'

# regular expressions
RE_EC = r'([0-9]\.{1}[^ ,;]*)'  # look for EC number format
RE_XTRA_SP = r' {2,}'  # look for 2+ spaces
RE_REC_VAL = r':(.*);'  # record values are between : and ; for specific records
RE_WORDS = r'\s*(\S[\S| ]*\S)\s*'

# directory/path values
SEP = os.sep
path_cwd = os.getcwd()
out_dir = path_cwd + SEP + 'pdb-output' + SEP
image_dir = path_cwd + SEP + 'images' + SEP
formatted_out = out_dir + 'format_out_all.txt'
formatted_out_small = out_dir + 'format_out_partial.txt'
formatted_basic = out_dir + 'basic_info_formatted.txt'
pdb_dir = out_dir + 'pdb_files' + SEP

# default values
DS = ' '  # default string
DF = 0.0  # default float
DI = 0  # default int
DL = []  # default list

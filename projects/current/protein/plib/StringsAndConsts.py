import os

CHUNK_SIZE = 6
PART_URL = "https://files.rcsb.org/view/"

# PDB file keys
K_CMP = 'COMPND'
K_EC = 'EC:'
K_EX_SYS = 'EXPRESSION_SYSTEM'
K_EX_TAX = 'EXPRESSION_SYSTEM_TAXID'
K_HEAD = 'HEADER'
K_ORG = 'ORGAN'
K_ORG_CMN = 'ORGANISM_COMMON'
K_ORG_SCI = 'ORGANISM_SCIENTIFIC'
K_ORG_TAX = 'ORGANISM_TAXID'
K_SRC = 'SOURCE'
K_TTL = 'TITLE'
K_ATM = 'ATOM'
K_HAT = 'HETATM'

# regular expressions
RE_EC = '([0-9]\.{1}[^ ,;]*)'  # look for EC number format
RE_XTRA_SP = ' {2,}'  # look for 2+ spaces
RE_REC_VAL = '(:.*;)'  # record values are between : and ; for specific records

# directory/path values
SEP = os.sep
cwd = os.getcwd()
out_dir = cwd + SEP + 'pdb-output' + SEP
formatted_out = out_dir + 'format_out_all.txt'
formatted_out_small = out_dir + 'format_out_partial.txt'
pdb_dir = out_dir + 'pdb_files' + SEP
sasa = '..' + SEP + '..' + SEP + 'dr_sasa_n-0.4b' + SEP + 'build' + SEP + 'dr_sasa'
sasa_dir = cwd + SEP + 'sasa-output' + SEP

try:
    from util import *
except ModuleNotFoundError:
    from sharedlib.util import *
    from sharedlib.regexlines import *

# protein code stuff
# ligand_info = get_json_data(FN_LIGAND_INFO, JKEY)
# ligand_codes = get_json_data(FN_LIGAND_CODES, JKEY)
# pdb_id_list = get_json_data(FN_PDB_IDS, JKEY)
# pdb_id_list_short = get_json_data(FN_PDB_IDS_SHORT, JKEY)
# pdb_id_single = get_json_data(FN_PDB_SINGLE, JKEY)
# protein_info = get_json_data(FN_PROTEIN, JKEY)
# struct_rings = get_json_data(FN_RINGS, JKEY)

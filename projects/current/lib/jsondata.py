try:
    from util import *
    from pathstrings import *
except ModuleNotFoundError:
    from lib.util import *
    from lib.pathstrings import *

flav_list = get_json_data(FN_FLAVS, JKEY)
flav_synonyms = get_json_data(FN_FLAV_SYN)
flav_relatives = get_json_data(FN_FLAV_REL)

# fill lists and dictionaries from json files
test_single = get_json_data(FN_SINGLE, JKEY)
test_short = get_json_data(FN_SHORT, JKEY)
test_med = get_json_data(FN_MED, JKEY)

# species codes for plants of interest
plant_full_list = get_json_data(FN_PLANT_LIST, JKEY)
plant_list = plant_full_list

# Dictionary/list of key-value pairs that defines the appropriate genus species for each code
plant_dict_reg = get_json_data(FN_PLANT_DICT, JKEY)
plant_dict_common = get_json_data(FN_PLANT_DICT_COMMON, JKEY)
plant_dict = plant_dict_reg

# list of pathways of interest
path_map_list = get_json_data(FN_PATH_LIST, JKEY)

# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
path_map_dict = get_json_data(FN_PATH_DICT, JKEY)

# protein code stuff
ligand_info = get_json_data(FN_LIGAND, JKEY)
ligand_codes = get_json_data(FN_LIGAND_CODES, JKEY)
pdb_id_list = get_json_data(FN_PDB_IDS, JKEY)
pdb_id_list_short = get_json_data(FN_PDB_IDS_SHORT, JKEY)
pdb_id_single = get_json_data(FN_PDB_SINGLE, JKEY)
protein_info = get_json_data(FN_PROTEIN, JKEY)
struct_planes = get_json_data(FN_PLANES)

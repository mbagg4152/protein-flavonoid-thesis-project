from lib.util import *
from lib.file_path_names import *

# fill lists and dictionaries from json files
test_single = get_json_data(json_single, j_key)
test_short = get_json_data(json_short, j_key)
test_med = get_json_data(json_med, j_key)

# main lists for the program
# species codes for plants of interest,  "mtr" cut due to code errors
full_list = get_json_data(f_species_list_json, j_key)

# Dictionary/list of key-value pairs that defines the appropriate genus species for each code
species_pairs_original = get_json_data(f_species_dict_json, j_key)
species_pairs_common = get_json_data(json_dict_common, j_key)
species_pairs = species_pairs_original
# list of pathways of interest
pathway_list = get_json_data(path_list_json, j_key)

# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
pathway_pairs = get_json_data(path_dict_json, j_key)

from lib.util import *
from lib.file_path_names import *

# fill lists and dictionaries from json files
test_single = get_json_data(fn_test_single, j_key)
test_short = get_json_data(fn_test_short, j_key)
test_med = get_json_data(fn_test_med, j_key)

# main lists for the program
# species codes for plants of interest,  "mtr" cut due to code errors
alpha_bact_list = get_json_data(fn_alpha_bact_list, j_key)
cyan_bact_list = get_json_data(fn_cyan_bact_list, j_key)
fungi_list = get_json_data(fn_fungi_list, j_key)
planc_bact_list = get_json_data(fn_planc_bact_list, j_key)
plant_list_new = get_json_data(fn_plant_list_new, j_key)
plant_list_orig = get_json_data(fn_plant_list, j_key)
protist_list = get_json_data(fn_protist_list, j_key)

full_list = plant_list_orig

# Dictionary/list of key-value pairs that defines the appropriate genus species for each code
alpha_bact_pairs = get_json_data(fn_alpha_bact_pairs, j_key)
cyan_bact_pairs = get_json_data(fn_cyan_bact_pairs, j_key)
fungi_pairs = get_json_data(fn_fungi_pairs, j_key)
planc_bact_pairs = get_json_data(fn_planc_bact_pairs, j_key)
plant_pairs_new = get_json_data(fn_plant_pairs_new, j_key)
plant_pairs_orig = get_json_data(fn_plant_pairs, j_key)
plant_pairs_plus = get_json_data(fn_plants_pairs_plus, j_key)
protist_pairs = get_json_data(fn_protist_pairs, j_key)

species_pairs = plant_pairs_orig

# list of pathways of interest
pathway_list = get_json_data(fn_path_list, j_key)

# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
pathway_pairs = get_json_data(fn_path_pairs, j_key)

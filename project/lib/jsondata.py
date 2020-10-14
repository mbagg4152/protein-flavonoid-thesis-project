from util import *
from pathstrings import *

flav_list = get_json_data(fn_flavs, j_key)

# fill lists and dictionaries from json files
test_single = get_json_data(fn_test_single, j_key)
test_short = get_json_data(fn_test_short, j_key)
test_med = get_json_data(fn_test_med, j_key)

# main lists for the program
# species codes for plants of interest,  "mtr" cut due to code errors
plant_list_orig = get_json_data(fn_plant_list, j_key)

full_list = plant_list_orig

# Dictionary/list of key-value pairs that defines the appropriate genus species for each code

plant_pairs = get_json_data(fn_plant_pairs, j_key)
plant_pairs_plus = get_json_data(fn_plants_pairs_plus, j_key)

species_pairs = plant_pairs

# list of pathways of interest
pathway_list = get_json_data(fn_path_list, j_key)

# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
pathway_pairs = get_json_data(fn_path_pairs, j_key)

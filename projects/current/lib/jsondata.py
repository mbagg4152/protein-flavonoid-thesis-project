from lib.util import *
from lib.pathstrings import *

flav_list = get_json_data(FN_FLAVS, JKEY)

# fill lists and dictionaries from json files
test_single = get_json_data(FN_TEST_SINGLE, JKEY)
test_short = get_json_data(FN_TEST_SHORT, JKEY)
test_med = get_json_data(FN_TEST_MED, JKEY)

# species codes for plants of interest
full_list = get_json_data(FN_PLANT_LIST, JKEY)
plant_list = test_short

# Dictionary/list of key-value pairs that defines the appropriate genus species for each code
plant_dict_reg = get_json_data(FN_PLANT_DICT, JKEY)
plant_dict_common = get_json_data(FN_PLANT_DICT_COMMON, JKEY)
plant_dict = plant_dict_reg

# list of pathways of interest
path_list = get_json_data(FN_PATH_LIST, JKEY)

# Dictionary/list of key-value pairs defining each pathway by each chemical they're responsible for
path_dict = get_json_data(FN_PATH_DICT, JKEY)

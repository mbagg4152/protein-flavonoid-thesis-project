import os
from parseutil import *

SEP = os.sep

flav_names = 'njson' + SEP + 'flav_names.json'
flav_rel = 'njson' + SEP + 'flav_related.json'
flav_syns = 'njson' + SEP + 'flav_syns.json'
ncomp_codes = 'njson' + SEP + 'chem_codes.json'
ncomp_dict = 'njson' + SEP + 'chem_codes_names.json'
nplant_codes = 'njson' + SEP + 'organism_codes.json'
nplant_dict = 'njson' + SEP + 'organism_codes_names.json'
plant_names_codes = 'njson' + SEP + 'plant_names_codes.json'

flav_list = get_json_data(flav_names)
flav_relatives = get_json_data(flav_rel)
flav_synonyms = get_json_data(flav_syns)
plant_dict = get_json_data(plant_names_codes)

cmp_codes = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_codes.txt'
comp_dict = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_dict.txt'
comp_in = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_info.txt'
org_codes = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'plant_codes.txt'
org_in = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'npass_species.txt'
org_out = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'parsed_species.txt'
pairs_in = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'plants_compounds.txt'
selected = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'selected.txt'

re_comp_id = r'(NPC\S*)'
re_comp_name = r'NPC\S*\s*([^\t\r\n]*)'
re_has_cas = r'<tr><td class="d1"><a href=information\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td><td ' \
             r'class="d1">([0-9]*-[0-9]*-[0-9]*)<\/td><td class="d1">([^<]*)<\/td><td class="d1"> '
re_no_cas = r'<tr><td class="d1"><a h.*n\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td><td ' \
            r'class="d1"><\/td><td class="d1">([^<]*)<\/td><td class="d1"> '
re_org_comp_id = r'-(NPC\S*)'
re_org_id = r'(NPO\S*)'
re_org_name = r'NPO[^ \t\r]*\s*(\S*\s*\S*)\s*\bSpecies'
re_org_name_pairs = r'(NPO\S*)\b-'


URL_KNAP = "'http://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word="  # knapsack partial URL

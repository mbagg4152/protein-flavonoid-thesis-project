import os

slash = os.sep  # get the right slash. / for linux, \ for windows

# output file names
file_apig = 'apigenin.txt'
file_bute = 'butein.txt'
file_cate = 'catechin.txt'
file_cyan = 'cyanidin.txt'
file_ecat = 'epicatechin.txt'
file_epig = 'epigallocatechin.txt'
file_erio = 'eriodictyol.txt'
file_gall = 'gallocatechin.txt'
file_geni = 'genistein.txt'
file_kaem = 'kaempferol.txt'
file_lute = 'luteolin.txt'
file_myri = 'myricetin.txt'
file_nari = 'naringenin.txt'
file_quer = 'quercetin.txt'
file_ec1 = 'EC-2-4-1-74.txt'
file_ec2 = 'EC-2-3-1-70.txt'
file_ec3 = 'EC-2-3-1-30.txt'
file_ec4 = 'EC-2-4-1-136.txt'
readme = slash + 'ReadMe.txt'

# output directories
chem_dir = slash + 'Chemical_Data'
fasta_dir = slash + 'FASTA_Data'
gene_dir = slash + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
info_dir = 'json-data' + slash
f_species_dict_json = info_dir + 'plant_names_codes.json'
f_species_list_json = info_dir + 'plant_codes.json'
j_key = 'obj'
path_dict_json = info_dir + 'path_names_codes.json'
path_list_json = info_dir + 'path_codes.json'
json_single = info_dir + 'test_single.json'
json_short = info_dir + 'test_short.json'
json_med = info_dir + 'test_med.json'
json_dict_common = info_dir + 'plant_names_codes_more.json'
json_addtl_plant_dict = info_dir + 'addtl_plant_names_codes.json'
json_addtl_plant_list= info_dir + 'addtl_plant_codes.json'

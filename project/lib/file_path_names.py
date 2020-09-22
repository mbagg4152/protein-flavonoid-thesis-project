import os

slash = os.sep  # get the right slash. / for linux, \ for windows

# output file names
fn_apig = 'apigenin.txt'
fn_bute = 'butein.txt'
fn_cate = 'catechin.txt'
fn_cyan = 'cyanidin.txt'
fn_ecat = 'epicatechin.txt'
fn_epig = 'epigallocatechin.txt'
fn_erio = 'eriodictyol.txt'
fn_gall = 'gallocatechin.txt'
fn_geni = 'genistein.txt'
fn_kaem = 'kaempferol.txt'
fn_lute = 'luteolin.txt'
fn_myri = 'myricetin.txt'
fn_nari = 'naringenin.txt'
fn_quer = 'quercetin.txt'
fn_ec1 = 'EC-2-4-1-74.txt'
fn_ec2 = 'EC-2-3-1-70.txt'
fn_ec3 = 'EC-2-3-1-30.txt'
fn_ec4 = 'EC-2-4-1-136.txt'
fn_readme = slash + 'ReadMe.txt'

# output directories
chem_dir = slash + 'Chemical_Data'
fasta_dir = slash + 'FASTA_Data'
gene_dir = slash + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
info_dir = 'json-data' + slash
j_key = 'obj'
fn_alpha_bact_list = info_dir + 'ap_bacteria_list.json'
fn_alpha_bact_pairs = info_dir + 'ap_bacteria_pairs.json'
fn_cyan_bact_list = info_dir + 'cyan_bact_list.json'
fn_cyan_bact_pairs = info_dir + 'cyan_bact_pairs.json'
fn_fungi_list = info_dir + 'fungi_codes.json'
fn_fungi_pairs = info_dir + 'fungi_names_codes.json'
fn_path_list = info_dir + 'path_codes.json'
fn_path_pairs = info_dir + 'path_names_codes.json'
fn_planc_bact_list = info_dir + 'planc_bacteria_list.json'
fn_planc_bact_pairs = info_dir + 'planc_bacteria_pairs.json'
fn_plant_list = info_dir + 'plant_codes.json'
fn_plant_list_new = info_dir + 'addtl_plant_codes.json'
fn_plant_pairs = info_dir + 'plant_names_codes.json'
fn_plant_pairs_new = info_dir + 'addtl_plant_names_codes.json'
fn_plants_pairs_plus = info_dir + 'plant_names_codes_more.json'
fn_protist_list = info_dir + 'protist_codes.json'
fn_protist_pairs = info_dir + 'protist_names_codes.json'
fn_test_med = info_dir + 'test_med.json'
fn_test_short = info_dir + 'test_short.json'
fn_test_single = info_dir + 'test_single.json'

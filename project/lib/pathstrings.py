import os

slash = os.sep  # get the right slash. / for linux, \ for windows

# output file names
fn_apig = 'apigenin.csv'
fn_bute = 'butein.csv'
fn_cate = 'catechin.csv'
fn_cyan = 'cyanidin.csv'
fn_ecat = 'epicatechin.csv'
fn_epig = 'epigallocatechin.csv'
fn_erio = 'eriodictyol.csv'
fn_gall = 'gallocatechin.csv'
fn_geni = 'genistein.csv'
fn_kaem = 'kaempferol.csv'
fn_isol = 'isoliquiritigenin.csv'
fn_lute = 'luteolin.csv'
fn_myri = 'myricetin.csv'
fn_nari = 'naringenin.csv'
fn_quer = 'quercetin.csv'
fn_ec23 = 'EC-2-4-1-74.csv'
fn_ec24 = 'EC-2-3-1-70.csv'
fn_ec25 = 'EC-2-3-1-30.csv'
fn_ec26 = 'EC-2-4-1-136.csv'
fn_readme = slash + 'ReadMe.txt'

# output directories
chem_dir = slash + 'Chemical_Data'
fasta_dir = slash + 'FASTA_Data'
gene_dir = slash + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
info_dir = 'json-data' + slash
j_key = 'obj'
fn_path_list = info_dir + 'path_codes.json'
fn_path_pairs = info_dir + 'path_names_codes.json'
fn_plant_list = info_dir + 'plant_codes.json'
fn_plant_pairs = info_dir + 'plant_names_codes.json'
fn_plants_pairs_plus = info_dir + 'plant_names_codes_more.json'
fn_test_med = info_dir + 'test_med.json'
fn_test_short = info_dir + 'test_short.json'
fn_test_single = info_dir + 'test_single.json'

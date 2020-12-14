from bioservices.kegg import KEGG
import datetime
import json
import os
import re
import sys
import threading
import urllib.error
import urllib.parse
import urllib.request


########################################################################################################################
# CONSTANTS
########################################################################################################################
SEP = os.sep  # get the right slash. / for linux, \ for windows

# output file names
FN_AGI, FN_BUN, FN_DEC, FN_EC = 'apigenin.csv', 'butein.csv', 'EC-2-3-1-70.csv', 'epicatechin.csv'
FN_EGT, FN_ERD, FN_GC, FN_GEN = 'epigallocatechin.csv', 'eriodictyol.csv', 'gallocatechin.csv', 'genistein.csv'
FN_GGT, FN_HCC, FN_HWB, FN_KMP = 'EC-2-4-1-74.csv', 'isoliquiritigenin.csv', 'cyanidin.csv', 'kaempferol.csv'
FN_KXN, FN_LU2, FN_MYC, FN_NAR = 'catechin.csv', 'luteolin.csv', 'myricetin.csv', 'naringenin.csv'
FN_QUER, FN_SOA, FN_V1G, FN_README = 'quercetin.csv', 'EC-2-3-1-30.csv', 'EC-2-4-1-136.csv', SEP + 'ReadMe.txt'

# output directories
CHEM_DIR, FASTA_DIR, GENE_DIR = SEP + 'Chemical_Data', SEP + 'FASTA_Data', SEP + 'Gene_Data'

# pathways for the json data  (and the name of the json object used in all files)
JSON_DIR = '..' + SEP + 'json_data' + SEP
FLAV_JSON, PROT_JSON, JKEY = JSON_DIR + 'flavonoid' + SEP, JSON_DIR + 'protein' + SEP, 'obj'
FN_FLAVS, FN_PATH_DICT = FLAV_JSON + 'flav_names.json', FLAV_JSON + 'flav_path_names_codes.json'
FN_PATH_LIST, FN_PLANT_DICT = FLAV_JSON + 'flav_path_codes.json', FLAV_JSON + 'flav_plant_names_codes.json'
FN_PLANT_DICT_COMMON, FN_PLANT_LIST = FLAV_JSON + 'flav_plant_names_codes_more.json', FLAV_JSON + 'flav_plant_codes.json'
FN_MED, FN_SHORT, FN_SINGLE = FLAV_JSON + 'flav_test_med.json', FLAV_JSON + 'flav_test_short.json', FLAV_JSON + 'flav_test_single.json'

# misc strings
NL, SP, NIX, CSV, GDATA = '\n', ' ', '', '.csv', 'Gene_data_'
URL_DBGET = 'https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+'

# keys for accessing dictionaries
O_KEY, E_KEY, N_KEY, P_KEY, G_KEY = 'ORTHOLOGY', 'EC', 'NTSEQ', 'PLANT', 'GENE'

# regular expressions
RE_NT_HEAD, RE_NT_SEQ, RE_SQ_BRACKETS = r'(&gt;.*)', r'([acgt]{10,})', r"(\[.*\])"
RE_ALPH, RE_EC, RE_KO = r'[^a-zA-Z]+', r'([ |:][0-9]*\.[0-9|-]*\.[0-9|-]*\.[0-9|-]*)', r'(\[KO:K[^\]]*\])'

# labels for each of the compounds
AGI, BUN, CAQ, EC, EGT = 'Apigenin', 'Butein', 'Catechol', 'Epicatechin', 'Epigallocatechin'
ERD, GC, GEN, HCC, HWB = 'Eriodictyol', 'Gallocatechin', 'Genistein', 'Isoliquiritigenin', 'Cyanidin'
KMP, KXN, LU2, MYC, NAR = 'Kaempferol', 'Catechin', 'Luteolin', 'Myricetin', 'Naringenin'
PYG, QUE, RCO, STL = 'Pyrogallol', 'Quercetin', 'Resorcinol', 'Resveratrol'
E01, E02, E03, E04, E05 = 'EC:4.3.1.24', 'EC:4.3.1.25', 'EC:1.14.14.91', 'EC:6.2.1.12', 'EC:2.3.1.170'
E06, E07, E08, E09, E10 = 'EC:2.3.1.133', 'EC:1.14.14.96', 'EC:1.14.13.-', 'EC:2.3.1.74', 'EC:5.5.1.6'
E11, E12, E13, E14, E15 = 'EC:1.14.20.5', 'EC:1.14.19.76', 'EC:1.14.14.81', 'EC:1.14.14.82', 'EC:1.14.11.9'
E16, E17, E17_2, E18, E19 = 'EC:1.14.20.6', 'EC:1.1.1.219', 'EC:1.1.1.219 1.1.1.234', 'EC:1.14.20.4', 'EC:1.3.1.77'
E20, E21, E22, GGT, DEC = 'EC:1.17.1.3', 'EC:1.14.14.87', 'EC:4.2.1.105', 'EC:2.4.1.74', 'EC:2.3.1.70'
SOA, V1G = 'EC:2.3.1.30', 'EC:2.4.1.136'

########################################################################################################################
# GLOBAL VARIABLES
########################################################################################################################
init_time = datetime.datetime.now()
kegg = KEGG()

list_all_genes = []  # will hold a list of all gene objects
list_all_plants = []  # will hold a list of all plant objects
list_fasta_ec = []  # will hold a list of the objects listing FASTA by EC number
list_genes_by_path = []  # will hold a list of objects which are collections of genes for each plant pathway
list_plant_paths = []  # will hold a list of all plant pathways
list_all_plant_matrix = []  # will hold the list after calculating the count for each EC for each plant

lock_access_ec = threading.Lock()  # will allow only one thread to access the list of ec numbers at any time
lock_access_plant = threading.Lock()  # will allow only one thread to access the list of plants at any time
lock_add_gene = threading.Lock()  # will allow only one thread to access the list of genes at any time
lock_kegg_get = threading.Lock()  # will allow only one thread to get an entry from KEGG at any time

path_chem = ' '  # will hold path of directory that contains the prediction output
path_cwd = ' '  # current working directory at initial runtime
path_fasta = ' '  # will hold path of directory that contains the fasta output
path_gene = ' '  # will hold path of directory that contains the gene data output
path_main = ' '  # will hold path of parent directory of program
thread_lim = 5  # max number of threads to be used in accessing data

path_map_dict = {}  # will hold dictionaries of path map numbers and their names
path_map_list = []  # will hold list of path map numbers
plant_dict = {}  # will hold dictionary of the plant codes and corresponding names, assigned based on test case
plant_dict_common = {}  # will hold dictionary of plant codes and species names + common names
plant_dict_reg = {}  # will hold dictionary of plant codes and species names
plant_full_list = []  # will hold full list of kegg plant codes
plant_list = []  # will hold the list of KEGG codes used in the program. change in init_setup based on test case
test_med = []  # will hold medium sized list of KEGG codes, used for testing only
test_short = []  # will hold  short list of KEGG codes, used for testing only
test_single = []  # will hold list of 1 KEGG code, used for testing only

data_lists = []  # will hold list of flavonoid objects

########################################################################################################################
# MAIN KEGG CODE - PREDICTIONS, MATRIX & FASTA OUTPUTS
########################################################################################################################
def main():
    """
    This is the main function of the file which calls specific functions in order when running and displays the total run
    time at the end of the code execution.
    """
    init_setup()
    run_path_parse()
    prediction()
    run_fill_matrix()
    run_build_fasta()

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\nRun time: ' + str(total_time))

def init_setup():
    """
    This is the initial setup function for the program. If the user supplies a directory name in the command line args,
    then that name will be used when outputting the data. If no name is supplied, then the data is outputted to the
    directory named data. This function also initializes certain lists and dictionaries and creates the list of
    pathway and plant codes based on the KEGG codes that can be found in the JSON files. After making the list, it
    creates a list of plant objects that will be used throughout the program.
    """
    global path_cwd, path_main, path_chem, path_fasta, path_gene, list_plant_paths, list_all_plants, plant_full_list, \
        test_med, test_short, test_single, plant_list, plant_dict_reg, plant_dict_common, plant_dict, path_map_list, \
        path_map_dict, data_lists
    decision = ''
    if len(sys.argv) > 1:  # if length of args is greater than one that means user supplied arg other than program name
        decision = sys.argv[1]  # read in the desired output directory as a commandline argument
    else:
        print("No directory name supplied, defaulting to 'data'. Supply name using 'python3 keggv2.py dir_name'")
        decision = 'data'
    path_cwd = os.getcwd() + SEP  # get directory where program is being run
    path_main = path_cwd + decision  # make the main/data dir path
    path_chem = path_main + CHEM_DIR  # make the chem data path for prediction output
    path_fasta = path_main + FASTA_DIR  # make the fasta data path for the fetched fasta data
    path_gene = path_main + GENE_DIR  # make the gene data path for the fetched gene data

    init_dirs(path_main, path_gene, path_fasta, path_chem)  # (util.py) initialize directories if they don't exist

    # fill the lists and dictionaries from the JSON files
    test_single = get_json_data(FN_SINGLE, JKEY)
    test_short = get_json_data(FN_SHORT, JKEY)
    test_med = get_json_data(FN_MED, JKEY)
    plant_full_list = get_json_data(FN_PLANT_LIST, JKEY)
    plant_list = plant_full_list
    plant_dict_reg = get_json_data(FN_PLANT_DICT, JKEY)
    plant_dict_common = get_json_data(FN_PLANT_DICT_COMMON, JKEY)
    plant_dict = plant_dict_reg
    path_map_list = get_json_data(FN_PATH_LIST, JKEY)
    path_map_dict = get_json_data(FN_PATH_DICT, JKEY)

    # initialize the list of flavonoid objects
    data_lists = [ChemData(AGI, [], FN_AGI), ChemData(BUN, [], FN_BUN), ChemData(EC, [], FN_EC),
                  ChemData(HWB, [], FN_HWB), ChemData(EC, [], FN_EC), ChemData(EGT, [], FN_EGT),
                  ChemData(ERD, [], FN_ERD), ChemData(GC, [], FN_GC), ChemData(GEN, [], FN_GEN),
                  ChemData(KMP, [], FN_KMP), ChemData(LU2, [], FN_LU2), ChemData(MYC, [], FN_MYC),
                  ChemData(NAR, [], FN_NAR), ChemData(QUE, [], FN_QUER), ChemData(HCC, [], FN_HCC)]

    list_plant_paths = [i + j for i in plant_list for j in path_map_list]  # make plant and path combination list

    for key in plant_dict:  # build list of plant objects
        tmp_plant = Plant(code=key, name=plant_dict[key])  # object made for current plant
        if not tmp_plant.is_in(list_all_plants): list_all_plants.append(tmp_plant)  # add if not already present in list

def run_path_parse():
    """
    This function breaks the list of plant pathways into different lists in order for different data to be processed
    at the same time using multithreading. Once all of the threads have finished, then the list of genes by path will
    be looped through in order to create both the gene data output files for each pathway and for the master file
    that contains all of the gene information.
    """
    global list_plant_paths, list_all_plants, thread_lim, path_gene
    sub_lists = list_partition(list_plant_paths, thread_lim)  # split the list into parts
    threads = []  # will hold the created threads

    for sub_list in sub_lists:
        thread = threading.Thread(target=path_parse, args=(sub_list,))  # create thread & pass the sublist to path_parse
        thread.start()  # start running the thread
        threads.append(thread)  # add new thread to list of threads
    for thread in threads: thread.join()  # wait for all of the threads to finish running

    master_output = ''
    master_gene = path_main + SEP + 'MasterList.csv'
    for item in list_genes_by_path:
        tmp_file_path = path_gene + SEP + item.path + CSV
        tmp_output = ''  # will hold output from this gene
        for gene in item.genes:
            tmp_output += gene.simple() + '\n'  # use the simple function to get a formatted string for this pathway
            master_output += gene.simple() + '\n'  # get formatted string for the master file
        write_append(tmp_file_path, tmp_output)  # write the file for the pathway
    write_append(master_gene, master_output)  # write the master file

def path_parse(paths):
    """
    Given a pathway for a specific plant, the program then passes it to KEGG where it retrieves the appropriate entry.
    The EC number, KO number and each GENE entry are parsed from the data sent back by KEGG and then are appropriately
    saved by updating the plant objects.
    """
    for path in paths:
        global list_all_plants, list_genes_by_path
        print(path)
        with lock_kegg_get:
            raw = kegg.get(path)  # get KEGG entry for pathway
            gene_entry = kegg.parse(raw)  # parse kegg entry into dictionary for easier access
            entry_dict = gene_entry.get(G_KEY)  # get the data for the dictionary key GENE
        if entry_dict is not None:
            plant_code = ''.join(re.split(RE_ALPH, path))  # get only the letters from the path
            plant_name = plant_dict.get(plant_code)  # get the plant name by accessing the plant dictionary
            with lock_add_gene: list_genes_by_path.append(PathGene(path=path))  # add new path to the list
            for key in entry_dict:
                try:
                    # find EC number in the entry using regular expressions then remove square brackets
                    ec_num = re.findall(RE_EC, entry_dict[key])

                    for i in range(0, len(ec_num)):
                        item = ec_num[i].replace('[', '').replace(']', '').replace(':', '').replace(' ', '')
                        ec_num[i] = 'EC:' + item

                    # get the orthology ID using regular expressions then remove square brackets
                    orthology = mult_replace(quick_fetch(RE_KO, entry_dict[key]), [('[', ''), (']', '')])

                    # remove the EC & KO using regular expressions in order to get the compound name
                    name = re.sub(RE_KO, '', (re.sub(RE_EC, '', entry_dict[key])))

                    # create new gene object using the information from kegg
                    tmp_gene = Gene(gene_id=key, plant=plant_name, ec_nums=ec_num, ortho=orthology,
                                    compound=name, path=path, plant_code=plant_code)
                    with lock_add_gene:
                        for gene_path in list_genes_by_path:
                            if gene_path.path == path: gene_path.genes.append(tmp_gene)  # add to list of genes for path
                    with lock_access_plant:
                        for index, plant in enumerate(list_all_plants):
                            if plant.name == tmp_gene.plant:
                                tmp_plant = plant
                                # add to list of plants genes if not already present
                                if not tmp_gene.is_in(tmp_plant.genes): tmp_plant.genes.append(tmp_gene)
                                # add to list of all genes if not already present
                                if not tmp_gene.is_in(list_all_genes): list_all_genes.append(tmp_gene)
                                tmp_plant.ec_nums.extend(ec_num)  # add to the plants list of ec numbers (dupes okay)
                                list_all_plants[index] = plant  # update the list of plants with modified plant object
                except IndexError: pass  # couldn't find items using regular expression findall

def prediction():
    """
    This is the function that goes through each plant, looks at the list of EC numbers then applies a function in order
    to determine whether or not the plant has the required EC numbers needed to synthesize each compound.
    """
    global list_all_plants
    for plant in list_all_plants:
        for chem_data in data_lists:
            if flav_check(chem_data.label, plant.ec_nums):  # passes check, has all of the flavonoids
                chem_data.plants.append(plant.name)  # add plant to flavonoids list

    # create the prediction output files for each flavonoid
    for key in data_lists:
        save_file([key.plants], key.file_name, path_chem)
        item_count = len(key.plants)
        print(key.label + ' predicted in ' + str(item_count) + ' entries. ' +
              'Data saved in ' + path_chem + SEP + key.file_name + '.')

def run_fill_matrix():
    """
    This function creates and outputs a 'matrix' relating to each species and EC number by running the fill_matrix
    function on multiple threads. For each gene entry containing a specific EC number, the program will increase the
    counter and display it at the end next to the appropriate EC number.
    """
    global list_all_plants
    for plant in list_all_plants: fill_count_matrix(plant)  # update the matrix using each plant

    out = ''
    for plant in list_all_plant_matrix:
        out += plant.name + ': '  # create beginning of line for current plant
        for count in plant.ec_counts:
            out += '[' + str(count.number) + ' count: ' + str(count.count) + '] '  # add the EC numbers and counts
        out += '\n'
    write_append(path_main + SEP + 'MasterCount.csv', out)  # write string to the output file

def fill_count_matrix(plant):
    """
    This function builds the ec counts for each list.
    """
    global list_all_plants
    tmp_plant = plant
    for num in tmp_plant.ec_nums:
        if tmp_plant.has_ec_count(num): tmp_plant.incr_ec_count(num)  # increment count if number is already present
        else:
            # create new count object for current EC number
            tmp_count = EcCounts(number=num, count=1)
            plant.ec_counts.append(tmp_count)
    with lock_access_plant: list_all_plant_matrix.append(tmp_plant)  # update count matrix

def run_build_fasta():
    """
    This function uses multithreading and the information gathered from running run_path_parse in order to get the
    FASTA/DNA sequence for each of the gene entries that were found. As before, the program parses the list after
    all threads are done and then created a FASTA file for each EC number and created the Master FASTA file.
    """

    sub_lists = list_partition(list_all_genes, thread_lim)
    threads = []
    print('getting data for ' + str(len(list_all_genes)) + ' genes')
    for sub_list in sub_lists:
        thread = threading.Thread(target=build_fasta, args=(sub_list,))  # create the thread
        thread.start()  # start thread execution
        threads.append(thread)  # append it to the list of threads

    for thread in threads: thread.join()  # wait for all of the threads to be done before moving on

    print('Starting to gather data for FASTA files...')
    master_fasta = path_main + SEP + 'MasterFASTA.csv'
    master_output = ''
    for item in list_fasta_ec:
        tmp_file_path = path_fasta + SEP + item.ec_name.replace('.', '-').replace(':', '') + CSV
        tmp_output = ''
        for entry in item.ec_entries:
            tmp_output += entry.simple() + '\n'  # add to string to be printed into specific EC file
            master_output += entry.simple() + '\n'  # add to string to be printed into master FASTA
        write_append(tmp_file_path, tmp_output)  # write the file for current EC number
    write_append(master_fasta, master_output)  # write the master FASTA file
    print('Done making the FASTA files.')

def build_fasta(genes):
    """
    This function uses the plant code and gene id in order to find the matching FASTA sequence using the appropriate
    dbget url. The pages are saved into memory as HTML and are parsed in order to extract the important information from
    the web page. After parsing, the FASTA sequences are added to EcFastaCollection objects in order to maintain
    proper association when writing all of the sequences out to files.
    """
    global list_fasta_ec
    for gene in genes:
        # using the plant code and gene id create a string formatted as code:gene
        combined = gene.plant_code.strip() + ':' + gene.gene_id.replace('(RAP-DB) ', '').strip()
        db_url = URL_DBGET + combined  # append code-gene string to the end of the dbget incomplete URL
        try:
            # read the html from the dbget url
            with urllib.request.urlopen(db_url) as db_site: url_data = db_site.read().decode('utf-8')
        except urllib.error.HTTPError or urllib.error.URLError as e:  # error getting html data
            print('Something went wrong with error ' + e)
            continue

        # get the header of the FASTA entry using regular expressions, &gt; is the HTML representation of >
        fasta_header = ''.join(re.findall(RE_NT_HEAD, url_data)).replace('&gt;', '>')
        fasta_body = ''.join(re.findall(RE_NT_SEQ, url_data))  # get the DNA sequence body using regular expressions
        full_fasta_entry = fasta_header + '\n' + fasta_body  # create FASTA entry string
        # create new entry object
        tmp_entry = FastaEcEntry(gene=gene.gene_id, plant=plant_dict.get(gene.plant_code), dna=full_fasta_entry)
        with lock_access_ec:
            for g in gene.ec_nums:
                tmp_ec = EcFastaCollection(ec_num=g, ec_entries=[tmp_entry])
                list_fasta_ec.append(tmp_ec)

########################################################################################################################
# USED WITH LOGICAL FUNCTIONS
########################################################################################################################
def flav_check(label, ec_list):
    """
    determine which logical function to call based on the label passed in using the label function pairs
    """
    options = {
        AGI: agi, BUN: bun, KXN: kxn, HWB: hwb, EC: ec, EGT: egt, ERD: erd, GC: gc, GEN: gen, KMP: kmp,
        LU2: lu2, MYC: myc, NAR: nar, QUE: que, GGT: ggt, DEC: dec, SOA: soa, V1G: v1g, HCC: bun
    }
    res = ''
    try: res = options[label](ec_list)
    except KeyError: res = False
    return res

def or_in(items, *args):
    """
    returns true if at least 1 parameter is in the passed in EC list
    """
    found = 0
    for a in args:
        if a in items: found += 1
    if found > 0: return True
    else: return False

def and_in(items, *args):
    """
    returns true only if all parameters are in the passed in EC list
    """
    for a in args:
        if a not in items: return False
    return True

########################################################################################################################
# LOGICAL FUNCTIONS: returns true if the elements in the list meet the required conditions.
########################################################################################################################
def tca(e): return or_in(e, E01, E02)  # cinnamic acid
def hc4(e): return tca(e) and (E03 in e)  # p-coumaric acid
def nca(e): return tca(e) and (E04 in e)  # cinnamoyl-coa
def wca(e): return (nca(e) and (E03 in e)) or (hc4(e) and (E04 in e))  # p-coumaroyl-coa
def nch(e): return wca(e) and (E09 in e)  # naringenin chalcone
def nar(e): return nch and (E10 in e)  # naringenin
def agi(e): return nar(e) and or_in(e, E11, E12)  # apigenin
def lu2(e): return agi(e) and (or_in(e, E13, E14))  # luteolin
def fca(e): return and_in(e, E06, E07) or (E08 in e)  # caffeoyl-coa
def erd(e): return (nar(e) and or_in(e, E13, E14)) or (fca(e) and (E09 in e))  # eriodictyol
def dhk(e): return nar(e) and (E15 in e)  # dihydrokaempferol
def kmp(e): return dhk(e) and (E16 in e)  # kaempferol
def dhq(e): return (dhk(e) and or_in(e, E13, E14)) or (erd(e) and (E15 in e))  # dihydroquercetin
def lhwb(e): return dhq(e) and or_in(e, E17, E17_2)  # leucocyanidin
def que(e): return (kmp(e) and or_in(e, E13, E14)) or (dhq(e) and (E16 in e))  # quercetin
def kxn(e): return lhwb(e) and (E20 in e)  # catechin
def hwb(e): return lhwb(e) and (E18 in e)  # cyanidin
def ec(e): return hwb(e) and (E19 in e)  # epicatechin
def dhm(e): return (dhq(e) and (E13 in e)) or (erd(e) and and_in(e, E13, E15))  # dihydromyricetin
def myc(e): return (dhm(e) and (E16 in e)) or (que(e) and (E13 in e))  # myricetin
def ldlm(e): return dhm(e) and or_in(e, E17, E17_2)  # leucodelphinidin
def dlm(e): return ldlm(e) and (E18 in e)  # delphinidin
def gc(e): return ldlm(e) and (E20 in e)  # gallocatechin
def egt(e): return dlm(e) and (E19 in e)  # epigallocatechin
def gen(e): return nar(e) and and_in(e, E21, E22)  # genistein
def bun(e): return wca(e) and (or_in(e, E05, E09))  # butein
def ggt(e): return GGT in e  # Glycosaminoglycan galactosyltransferase
def dec(e): return DEC in e  # Deleted entry
def soa(e): return SOA in e  # Serine O-acetyltransferase
def v1g(e): return V1G in e  # vanillate 1-glucosyltransferase

########################################################################################################################
# UTILITIES
########################################################################################################################
def get_json_data(file_name, key):
    """
    This function uses the python JSON library in order to parse JSON files into usable python objects. Can return
    lists or dictionaries, depending on the JSON file's structure.
    """
    data = ''
    try:
        with open(file_name) as jsonFile: data = json.load(jsonFile)
    except FileNotFoundError:
        file_name = '..' + SEP + file_name
        with open(file_name) as jsonFile:
            data = json.load(jsonFile)
    return data[key]

def list_partition(seq, num):
    """
    This function takes in a list and then splits it into as many parts as specified by parameter num
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def save_file(lists_to_write, output_dir, current):
    """
    This function takes a list or list of lists and then writes its contents to a file.
    """
    os.chdir(current)

    writedoc = open(output_dir, 'w')  # open file to be written
    for line in lists_to_write:
        for item in line:
            item = str(item).replace(NL, NIX)  # removes the new lines in each list of list
            if item == NIX:   writedoc.write('-')  # if the list in the list of list is empty writes a dash
            else: writedoc.write(item)  # write the entry in the list of lists to the file
            writedoc.write(', ')  # tab delineated; use "," for csv files
    writedoc.write(NL)
    writedoc.close()

def write_append(path, content, write_over=None):
    """
    This function takes a file name and the contents to be written to a file. If the file doesn't exist, it is
    created then written to. If it does exist, then it is appended to. The optional arg write_over is used for when
    the file exists and the content needs to be written over.
    """
    try:
        file = open(path, 'x')
        file.close()
        file = open(path, 'w')
        file.write(content)
        file.close()
    except FileExistsError:
        if write_over:
            file = open(path, 'w')
            file.write(content)
            file.close()
        else:
            file = open(path, 'a')
            file.write(content)
            file.close()

def init_dirs(main_dir, gene, fasta, chem):
    """
    Initializes directories for keggv2.py
    """
    # replaced WindowsError with OSError for more general usage. try to make data directories and handle any errors
    try: os.mkdir(main_dir)
    except OSError: pass
    try: os.mkdir(gene)
    except OSError or FileExistsError: pass
    try: os.mkdir(fasta)
    except OSError or FileExistsError: pass
    try: os.mkdir(chem)
    except OSError or FileExistsError: pass

def quick_fetch(pattern, line):
    """
    Fetch one item from re.findall and return as a string.
    """
    out = ''
    try: out = re.findall(pattern, line)[0]
    except IndexError: out = ''
    return out

def mult_replace(line, pairs):
    """
    This function makes multiple string replacements.
    """
    for pair in pairs: line = line.replace(pair[0], pair[1])
    return line

########################################################################################################################
# OBJECTS
########################################################################################################################
class ChemData:
    """
    This class holds the data for each flavonoid. The objects are initialized with their file name and label and only
    later in the program, their empty list of plants will be filled.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.label: string that contains the flavonoids name
    self.plants: list of plants predicted to produce the flavonoid
    self.file_name: string that holds the flavonoids output file name
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    """
    def __init__(self, label: str, plants: [str], file_name: str):
        self.plants = plants
        self.label = label
        self.file_name = file_name
    def __eq__(self, other):
        return self.plants == other.plants and self.label == other.label and self.file_name == other.file_name
    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False

class EcFastaCollection:
    """
    This object is used to hold the associated FASTA entries for any given EC number.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.ec_name: the EC number & name used when writing the file
    self.ec_entries: the list of associated FASTA entries (FastaEcEntry objects)
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    """
    def __init__(self, ec_num=None, ec_entries=None):
        self.ec_name = ec_num if ec_num is not None else ' '
        self.ec_entries = ec_entries if ec_entries is not None else []
    def __eq__(self, other): return self.ec_name == other.ec_name and self.ec_entries == other.ec_entries
    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False

class EcCounts:
    """
    This object is a property of the Plant object and is used to hold each EC number and the number of times it occurs
    in gene entries of a given plant.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.number: the EC number
    self.count: number of times that the EC number shows up in gene entries.
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    """
    def __init__(self, number=None, count=None):
        self.number = number if number is not None else ' '
        self.count = count if number is not None else 0

class FastaEcEntry:
    """
    This object is a property of EcFastaCollection and contains the information for a specific FASTA entry.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.gene_id: the gene ID associated with the sequence
    self.plant: the plant that the gene is from
    self.dna_seq: the dna sequence/FASTA entry for the specific gene
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    simple: returns a formatted string
    """
    def __init__(self, gene=None, dna=None, plant=None):
        self.gene_id = gene if gene is not None else ' '
        self.dna_seq = dna if dna is not None else ' '
        self.plant = plant if dna is not None else ' '
    def __eq__(self, other):
        return self.gene_id == other.gene_id and self.dna_seq == other.dna_seq and self.plant == other.plant
    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False
    def simple(self): return self.dna_seq

class Gene:
    """
    This object holds data gathered from KEGG for each plant's pathway (like aip00491).
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.gene_id: the ID of the gene from a plant
    self.plant: the scientific name of the plant that has this gene
    self.plant_code: the KEGG code for the plant
    self.compound: the compound name listed in the entry
    self.ec_nums: the list of EC numbers found in the entry
    self.ortho: the KEGG orthology code for the compound
    self.path: the pathway where the gene was found
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    simple: returns a formatted string that contains information from the object
    no_plant: same as simple, but without including the plant name
    """
    def __init__(self, gene_id=None, plant=None, compound=None, ec_nums=None, ortho=None, path=None, plant_code=None):
        self.gene_id = gene_id if gene_id is not None else ' '
        self.plant = plant if plant is not None else ' '
        self.plant_code = plant_code if plant_code is not None else ' '
        self.compound = compound if compound is not None else ' '
        self.ec_nums = ec_nums if ec_nums is not None else []
        self.ortho = ortho if ortho is not None else ' '
        self.path = path if path is not None else ' '
    def __eq__(self, other):
        return self.gene_id == other.gene_id and self.plant_code == other.plant_code and self.ec_nums == other.ec_nums
    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False
    def simple(self):
        return self.plant + ', ' + self.gene_id + ', ' + self.compound + ', ' + str(self.ec_nums) + ', ' + self.ortho
    def no_plant(self):
        return self.gene_id + ', ' + self.compound + ', ' + self.ec_nums + ', ' + self.ortho

class PathGene:
    """
    This object is used to hold Gene objects in a way such that they are sorted by the pathway from which they were
    found.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.path: the pathway that resulted in the gene entry
    self.genes: the list of gene entries from this pathway
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    """
    def __init__(self, path=None, genes=None):
        self.path = path if path is not None else ' '
        self.genes = genes if genes is not None else []
    def __eq__(self, other):
        return self.path == other.path and self.genes == other.genes
    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False

class Plant:
    """
    This object holds information about each plant used in the program. The plant objects are initialized with their
    scientific name and their code and then have different information added later.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.name: scientific name of the plant
    self.code: KEGG code for the plant
    self.genes: the gene entries for the plant
    self.ec_nums: the EC numbers parsed from the plants gene entries
    self.flavonoids: the list of flavonoids that the plant could potentially produce
    self.ec_counts: list of objects that hold the number of times each EC number appears
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    has_ec_count: used to determine whether or not a specific EC number is already in the list of EC counts
    incr_ec_count: used to increase the count for the EC count objects.
    """
    def __init__(self, name=None, code=None, genes=None, ec_nums=None, flavonoids=None, ec_counts=None):
        self.name = name if name is not None else ' '
        self.code = code if code is not None else ' '
        self.genes = genes if genes is not None else []
        self.ec_nums = ec_nums if ec_nums is not None else []
        self.flavonoids = flavonoids if flavonoids is not None else []
        self.ec_counts = ec_counts if ec_counts is not None else []
    def __eq__(self, other):
        return self.name == other.name and self.code == other.code and self.genes == other.genes
    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False
    def simple(self):
        gstr = ''
        for gene in self.genes: gstr += gene.no_plant() + ' || '
        return self.name + ', ' + self.code + ', ' + gstr + ', ' + str(self.ec_nums) + ', ' + str(self.flavonoids)
    def has_ec_count(self, ec_number):
        for enum in self.ec_counts:
            if enum.number == ec_number: return True
        return False
    def incr_ec_count(self, ec_number):
        for enum in self.ec_counts:
            if enum.number == ec_number: enum.count += 1

if __name__ == '__main__':
    main()

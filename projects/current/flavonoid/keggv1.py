import datetime
import os
import re
import sys
import threading
import time
import warnings
import logging.config

sys.path.append(os.getcwd().replace(os.sep + 'flavonoid', ''))

from bioservices.kegg import KEGG
from lib.jsondata import *
from lib.datatypes import *
from lib.pathstrings import *
from lib.compoundinfo import *

init_time = datetime.datetime.now()
kegg = KEGG()

chem_path = ''
path_cwd = ''
fasta_path = ''
gene_path = ''
path_main = ''

count_matrix = [[]]
dna_dict = {}
enz_class_list = []
fasta_by_enz_class = []
fasta_by_class = {}
sem_master = threading.Semaphore()
sem_parsing = threading.Semaphore()
sem_write = threading.Semaphore()
sem_dna = threading.Semaphore()
lock_master = threading.Lock()
lock_parsing = threading.Lock()
lock_write = threading.Lock()
lock_dna = threading.Lock()
locks = []
master_gene_list = []
master_list = []
master_uniq = []
path_plant_list = []
plant_flavs = []
species_list = plant_list
thread_parsing_data = []
tmp_data_holder = []


def main():
    global path_plant_list
    path_plant_list = [i + j for i in species_list for j in path_list]
    start()
    master_pathway_parser()
    make_matrix_and_counts()
    make_fasta()
    write_readme(path_main, FN_README, init_time, fasta_path, gene_path)
    finish_up()


def start():
    global chem_path, fasta_path, gene_path, path_main, path_cwd
    decision = ''
    if len(sys.argv) > 1: decision = sys.argv[1]
    else:
        print("No directory name supplied in args, defaulting to directory 'data'. "
              "Supply directory name in terminal using 'python3 keggv1.py dir_name'")
        decision = 'data'

    # create sub dirs for gene, FASTA & chemical data
    cwd = os.getcwd() + SEP
    main_dir = cwd + decision
    chem_path = main_dir + CHEM_DIR
    fasta_path = main_dir + FASTA_DIR
    gene_path = main_dir + GENE_DIR

    # replaced WindowsError with OSError for more general usage. try to make data directories and handle any errors
    try: os.mkdir(main_dir)
    except OSError: pass
    try: os.mkdir(gene_path)
    except OSError or FileExistsError: pass
    try: os.mkdir(fasta_path)
    except OSError or FileExistsError: pass
    try: os.mkdir(chem_path)
    except OSError or FileExistsError: pass


def master_pathway_parser():
    global master_list, master_uniq, locks
    print('- parsing list of species & pathways...')
    threads = []
    chunked = list(list_partition(path_plant_list, CHUNK_SIZE))  # chunk up list to be run on multiple threads
    for i in range(0, len(chunked)):  # create semaphores and fill list of lists

        thread_parsing_data.append([])
        tmp_data_holder.append([])

    chunk_index = 0
    for chunks in chunked:  # run parse helper using the chunks on different threads
        t = threading.Thread(target=parse_helper, args=(chunks, chunk_index))
        t.start()
        threads.append(t)
        chunk_index += 1

    for t in threads: t.join()  # wait for threads to finish their tasks

    master_uniq = remove_dupes(master_list)
    master_uniq = list(filter(None, master_uniq))

    # removes false values and iterates through the list of lists
    for count in range(0, len(master_uniq)):  master_uniq[count] = list(filter(None, master_uniq[count]))


def parse_helper(plant_paths, sem_index):
    global master_list, tmp_data_holder, thread_parsing_data
    for path in plant_paths:
        with lock_parsing: tmp_data_holder[sem_index] = gene_pathway_data(path)

        with lock_write:
            thread_parsing_data[sem_index].extend(tmp_data_holder[sem_index])
            try: save_file(tmp_data_holder[sem_index], GDATA + path + CSV, gene_path)
            except AttributeError: return

    with lock_master: master_list.extend(thread_parsing_data[sem_index])


# function that fetches the required data
def gene_pathway_data(pathway_id):
    print(pathway_id)
    raw = kegg.get(pathway_id)
    gd = kegg.parse(raw)
    # print('data from gene_pathway_data:' + str(raw))
    line_count = 0
    gene_lines = []
    fetched_genes = gd.get('GENE')
    if fetched_genes is not None:
        gene_vals = fetched_genes.values()
        for gv in gene_vals: gene_lines.append(gv)
        for gene in gene_lines:  # this section makes a list of lists that are appreciatively separated
            split_sig = '^*^'
            gene = gene.replace('  ', split_sig).replace(';', split_sig).replace('[', split_sig).replace(']', '')
            gene = gene.split(split_sig)
            alpha_only = NIX
            for char in pathway_id:
                if char.isalpha(): alpha_only += char
            gene.insert(0, plant_dict[alpha_only])
            g_count = 0
            for g in gene:  # cleans up list of lists of extra spaces at the beginning and end of each list
                gene[g_count] = g.strip()
                g_count += 1
            gene_lines[line_count] = gene
            line_count += 1  # iterates through each list in the entry
    return gene_lines


def make_matrix_and_counts():
    global master_gene_list, count_matrix

    ec_list = unique_element_list(master_uniq, 'last')
    count_matrix = [['Species']]  # creating the matrix and adding up the counts
    count_matrix[0].extend(ec_list)  # adds the Unique EC numbers to the end of the matrix

    for sp in species_list:  count_matrix.append([sp])  # first item in each row (but first) is the matrix

    cols = 0
    for outer in count_matrix[0]:
        if cols != 0:  # first column isn't actually an EC#,
            ecc = outer
            rows = 0
            for inner in count_matrix:  # for each species (using the species code)
                if rows != 0:  # first row isn't actually a species
                    species = plant_dict[inner[0]]
                    counter3 = 0
                    for unique in master_uniq:  # iterate over the culled master list to check for matching sets
                        if unique[0] == species and unique[len(unique) - 1] == ecc: counter3 += 1
                    count_matrix[rows].append(str(counter3))
                rows += 1
        cols += 1

    # change master count to be actual species:
    count = 0
    for cm in count_matrix:
        if count != 0: count_matrix[count][0] = plant_dict[cm[0]]  # replaces species code with genus specie names
        count += 1

    # Make Master Files
    print('- making master files: no dupes & count matrix')
    save_file(master_uniq, 'Master_List.csv', path_main)
    save_file(count_matrix, 'Master_Count.csv', path_main)

    # make a master fasta file
    print('- about to make master FASTA')
    swapped_order = {v: k for k, v in plant_dict.items()}  # reverses dictionary keys and values
    master_gene_list = []

    # combines species codes and gene numbers in a list to be used for the master fasta function
    for uniq in master_uniq: master_gene_list.append(swapped_order[uniq[0]] + ':' + uniq[1])
    # print(str(master_gene_list))


def make_fasta():
    print('- saving master fasta...')
    master_fasta = get_master_fasta(master_gene_list)
    # print(str(master_fasta))
    # basic_write(fasta_path + SEP + 'Master_FASTA.csv', 'w', NIX)
    master_out = ''
    print('- looping through master fasta...')
    for key in master_fasta:
        tmp_dict = dna_dict[key]  # get dictionary for the current key (gene name)
        this_ec = tmp_dict[E_KEY]  # get ec number using the key EC
        # make the line to be added to the dictionary & then printed
        tmp_line = '>' + str(key) + ' ' + str(tmp_dict[P_KEY]) + ' ' + str(this_ec) + '\n' + str(tmp_dict[N_KEY][0])
        master_out += tmp_line + '\n'
        this_ec = this_ec.replace('[', NIX).replace(']', NIX)

        ec_flag = False

        tmp_data = ''
        try:
            fasta_by_class.get(this_ec)
            fasta_by_class[this_ec].append(tmp_line)
        except KeyError: fasta_by_class[this_ec] = [tmp_line]

    write_append(fasta_path + SEP + 'Master_FASTA.csv', master_out)

    print('- creating fasta files by ec number...')
    for key in fasta_by_class:  # make FASTA files by EC number
        name = key.replace('.', '_').replace(SP, NIX).replace(':', NIX)
        out = ''
        for item in fasta_by_class[key]: out += item + '\n'
        # basic_write(fasta_path + SEP + name + CSV, 'w', out)
        try:
            file = open(fasta_path + SEP + name + CSV, 'x')
            file.close()
        except FileExistsError: pass
        file = open(fasta_path + SEP + name + CSV, 'a')
        file.write(out)
        file.close()


def get_master_fasta(gene):
    print('- fetching data for master FASTA...')
    global master_gene_list, dna_dict
    threads = []
    partitions = list(list_partition(master_gene_list, CHUNK_SIZE))
    for part in partitions:
        # print(str(part))
        try:
            thread = threading.Thread(target=fasta_helper, args=(part,))
            thread.start()
            threads.append(thread)
        except threading.ThreadError:
            print('experienced a thread err')
            continue

    for thread in threads: thread.join()

    # os.system('echo ' + str(dna_dict) + ' > dna-dict.txt')
    return dna_dict


def fasta_helper(gene_list):
    global dna_dict

    for gene in gene_list:  # gene takes form kegg_code:gene_name
        gene = gene.replace('(RAP-DB) ', NIX)
        raw = kegg.get(gene)  # get kegg gene info
        if is_http_error(raw):
            for i in range(0, RETRY):
                raw = kegg.get(gene)
                if is_http_error(raw) or raw is None or raw == '' or raw is int:
                    if i == RETRY - 1: print('Made attempts to access data but failed, skipping entry for ' + gene)
                    else:
                        print('Got an HTTP error code from KEGG ... taking a small nap')
                        time.sleep(1)
                else: continue
        raw_info = kegg.parse(raw)  # turn kegg gene data into dictionary
        ortho_line = raw_info.get(O_KEY)  # get data from dictionary using the key ORTHOLOGY
        split_gene = gene.split(':')
        org_code = split_gene[0].strip()  # kegg code is first half
        gene_name = split_gene[1].strip()  # gene name is second half
        organism_name = plant_dict.get(org_code, NIX)  # get plant's scientific name from the plant dictionary

        # no orthology data, skip to next list item
        if ortho_line is None or isinstance(ortho_line, str):
            print('>>>passing, is none or string')
            continue

        ortho_val = ''  # will hold orthology values
        # extra safeguard - if data is not a dictionary then it is incorrect --> skip to next gene
        try: ortho_val = ''.join(ortho_line.values())
        except AttributeError:
            print('>>>passing, attribute err')
            continue

        parsed_ec = re.findall(RE_SQ_BRACKETS, ortho_val)  # use regular expression to find EC number, which is in []
        if len(parsed_ec) == 0 or len(parsed_ec[0]) < 1:
            print('>>>passing, no ec')
            continue  # no EC number, skip to next gene
        ec_num = parsed_ec[0]  # regular expression findall returns a list, expecting only one element
        ntseq = raw_info.get('NTSEQ', NIX)  # get values from kegg dictionary using NTSEQ as the key
        joined_dna_seq = [ntseq.replace(SP, NIX)]  # remove spaces from sequence

        with lock_dna:
            tmp_entry = {E_KEY: ec_num, N_KEY: joined_dna_seq, P_KEY: organism_name}
            dna_dict[gene_name] = tmp_entry
            os.system('echo ' + str(dna_dict) + ' >> dna-dict.txt')









def finish_up():
    master_ec_list = [['species', 'EC#s']]
    counter = 0
    print('- filling master matrix...')
    for entry in count_matrix:  # for each species
        species_ec_list = []
        if counter == 0: pass
        else:
            counter2 = 0
            for j in entry:  # for EC in species
                if counter2 == 0: species_ec_list.append(j)
                else:
                    if str(j) == "0": pass
                    else: species_ec_list.append(count_matrix[0][counter2])
                counter2 += 1
        master_ec_list.append(species_ec_list)
        counter += 1
    master_ec_list = master_ec_list[1:]

    print('- looping through master ec list...')
    for entry in master_ec_list:
        if len(entry) > 0:
            # tmp_entry = Species(entry[0], 0, [])
            for chem_data in data_lists:
                if flav_check(chem_data.label, entry):
                    chem_data.plants.append([entry[0]])
                    # tmp_entry.flavonoids.append(chem_data.label)
                    # tmp_entry.count += 1
            # plant_flavs.append(tmp_entry)

    for key in data_lists:
        save_file(key.plants, key.file_name, chem_path)
        item_count = len(key.plants)
        print(key.label + ' predicted in ' + str(item_count) + ' entries. ' +
              'Data saved in ' + chem_path + SEP + key.file_name + '.')

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken: ' + str(total_time))


main()

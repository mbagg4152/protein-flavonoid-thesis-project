import bioservices
import datetime
import os
import sys
import threading
import urllib.error
import urllib.parse
import urllib.request

sys.path.append(os.getcwd().replace(os.sep + 'flavonoid', ''))
from bioservices.kegg import KEGG
from lib.compoundinfo import *
from lib.datatypes import *
from lib.jsondata import *
from lib.pathstrings import *

"""
Global Variables used throughout the program
- init_time records the time the program starts
- kegg is the variable used to access the KEGG functions
- the lists are used when processing and accessing data and are also used in making the files
- The locks are to be used in the sections of code that use multithreading. Global variables are much easier
  to use when using multithreading but they also need protection. Specific variables must only be accessed
  one at a time or else unforeseen issues may show up.
- The path variables will hold the paths for the output directories based on the command line args (or lack 
  thereof) and where the program is being run.
"""
init_time = datetime.datetime.now()
kegg = KEGG()
list_all_genes = []
list_all_plants = []
list_fasta_ec = []
list_genes_by_path = []
list_plant_paths = []
list_all_plant_matrix = []

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


def main():
    """
    This is the main function of the file which calls specific functions in order when running and displays the total run
    time at the end of the code execution.
    """
    init_setup()
    run_path_parse()
    run_build_fasta()
    run_fill_matrix()
    prediction()

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\nRun time: ' + str(total_time))


def init_setup():
    """
    This is the initial setup function for the program. If the user supplies a directory name in the command line args,
    then that name will be used when outputting the data. If no name is supplied, then the data is outputted to the
    directory named data. This function also creates the list of pathway and plant codes based on the KEGG codes that
    can be found in the JSON files. After making the list, it creates a list of plant objects that will be used
    throughout the program.
    """
    global path_cwd, path_main, path_chem, path_fasta, path_gene, list_plant_paths, list_all_plants
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

    list_plant_paths = [i + j for i in plant_list for j in path_list]  # make plant and path combination list

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

    total_out = ''
    master_gene = path_main + SEP + 'MasterList.csv'
    for gene_path in list_genes_by_path:
        tmp_file = path_gene + SEP + gene_path.path + CSV
        out = ''  # will hold output from this gene
        for gene in gene_path.genes:
            out += gene.simple() + '\n'  # use the simple function to get a formatted string for this pathway
            total_out += gene.simple() + '\n'  # get formatted string for the master file
        basic_write(tmp_file, 'w', out)  # write the file for the pathway
    basic_write(master_gene, 'w', total_out)  # write the master file


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
        t = threading.Thread(target=build_fasta, args=(sub_list,))  # create the thread
        t.start()  # start thread execution
        threads.append(t)  # append it to the list of threads

    for t in threads: t.join()  # wait for all of the threads to be done before moving on

    master_fasta = path_main + SEP + 'MasterFASTA.csv'
    fasta_out = ''
    for fasta_ec in list_fasta_ec:
        tmp_name = path_fasta + SEP + fasta_ec.ec_name.replace('.', '-').replace(':', '') + CSV
        out = ''
        for entry in fasta_ec.ec_entries:
            out += entry.simple() + '\n'  # add to string to be printed into specific EC file
            fasta_out += entry.simple() + '\n'  # add to string to be printed into master FASTA
        basic_write(tmp_name, 'w', out)  # write the file for current EC number
    basic_write(master_fasta, 'w', fasta_out)  # write the master FASTA file

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
            out += '[' + count.number + ' count: ' + str(count.count) + '] '  # add the EC numbers and counts
        out += '\n'
    basic_write(path_main + SEP + 'MasterCount.csv', 'w', out)  # write string to the output file

def path_parse(paths):
    for path in paths:
        global list_all_plants, list_genes_by_path
        print(path)
        with lock_kegg_get:
            raw = kegg.get(path)
            gene_entry = kegg.parse(raw)  # parse kegg entry into dictionary for easier access
            entry_dict = gene_entry.get(G_KEY)  # get the data for the dictionary key GENE
        if entry_dict is not None:
            plant_code = ''.join(re.split(RE_ALPH, path))  # get only the letters from the path
            plant_name = plant_dict.get(plant_code)  # get the plant name by accessing the plant dictionary
            with lock_add_gene: list_genes_by_path.append(PathGene(path=path))  # add new path to the list
            for key in entry_dict:
                try:
                    # find EC number in the entry using regular expressions then remove square brackets
                    ec_num = mult_replace(quick_fetch(RE_EC, entry_dict[key]), [('[', ''), (']', '')])

                    # pass if this EC number is not relevant to the predictions
                    if ec_num not in ec_nums_of_interest: continue

                    # get the orthology ID using regular expressions then remove square brackets
                    orthology = mult_replace(quick_fetch(RE_KO, entry_dict[key]), [('[', ''), (']', '')])

                    # remove the EC & KO using regular expressions in order to get the compound name
                    name = re.sub(RE_KO, '', (re.sub(RE_EC, '', entry_dict[key])))

                    # create new gene object using the information from kegg
                    tmp_gene = Gene(gene_id=key, plant=plant_name, ec_num=ec_num, k_ortho=orthology,
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
                                tmp_plant.ec_nums.append(ec_num)  # add to the plants list of ec numbers (dupes okay)
                                list_all_plants[index] = plant  # update the list of plants with modified plant object
                except IndexError: pass  # couldn't find items using regular expression findall

def build_fasta(genes):
    global list_fasta_ec
    for gene in genes:
        combined = gene.plant_code.strip() + ':' + gene.gene_id.replace('(RAP-DB) ', '').strip()  # code:geneid
        db_url = DBGET_URL + combined
        try:
            with urllib.request.urlopen(db_url) as u: url_data = u.read().decode('utf-8')
        except urllib.error.HTTPError or urllib.error.URLError as e:
            print('got err ' + e)
            continue
        fasta_header = ''.join(re.findall(RE_NT_HEAD, url_data)).replace('&gt;', '>')  # change name
        fasta_body = ''.join(re.findall(RE_NT_SEQ, url_data))
        full_fasta_entry = fasta_header + '\n' + fasta_body
        tmp_entry = FastaEcEntry(gene=gene.gene_id, plant=plant_dict.get(gene.plant_code), dna=full_fasta_entry)
        with lock_access_ec:
            count = 0
            for i in range(0, len(list_fasta_ec)):
                if list_fasta_ec[i].ec_name == gene.ec_num:
                    tmp = list_fasta_ec[i]
                    if not tmp_entry.is_in(tmp.ec_entries): tmp.ec_entries.append(tmp_entry)
                    count += 1
            if count == 0:
                tmp_ec = EcFastaCollection(ec_num=gene.ec_num, ec_entries=[tmp_entry])
                if not tmp_ec.is_in(list_fasta_ec): list_fasta_ec.append(tmp_ec)

def prediction():
    global list_all_plants
    for plant in list_all_plants:
        for chem_data in data_lists:
            if flav_check(chem_data.label, plant.ec_nums):
                chem_data.plants.append(plant.name)

    for key in data_lists:
        save_file([key.plants], key.file_name, path_chem)
        item_count = len(key.plants)
        print(key.label + ' predicted in ' + str(item_count) + ' entries. ' +
              'Data saved in ' + path_chem + SEP + key.file_name + '.')

def fill_count_matrix(plant):
    global list_all_plants
    tmp_plant = plant
    # with lock_access_plant: replace_at = list_all_plants.index(plant)
    for num in tmp_plant.ec_nums:
        replace_num = tmp_plant.ec_nums.index(num)
        if tmp_plant.has_ec_count(num):
            tmp_plant.incr_ec_count(num)
        else:
            tmp_count = EcCounts(number=num, count=1)
            plant.ec_counts.append(tmp_count)
    with lock_access_plant: list_all_plant_matrix.append(tmp_plant)

if __name__ == '__main__':
    main()

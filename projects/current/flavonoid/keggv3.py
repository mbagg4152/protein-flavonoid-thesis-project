from flib.fconstants import *
from flib.data_types import *
import flib.prediction_logic as predict
from urllib import request
from urllib.error import HTTPError, URLError
import datetime
import sys
import threading

sys.path.append(os.getcwd().replace(os.sep + 'flavonoid', ''))
from sharedlib.regexlines import *

try:
    import bioservices
    from bioservices.kegg import KEGG
except ImportError as e:
    print('Program needs bioservices in order to work. In the terminal, try:\n'
          '`pip3 install bioservices` or `pip install bioservices`')
    print(e)
    exit(1)

####################
# Global Variables #
####################
init_time = datetime.datetime.now()  # get time of program execution
kegg = KEGG()  # used to access KEGG's functions from bioservices

# The lists are used for processing/accessing data as well as making output files.
list_all_genes = []
list_all_plant_matrix = []
list_all_plants = []
list_fasta_ec = []
list_genes_by_path = []
list_plant_paths = []

# The locks are used to protect the values of the global variables when using multithreading.
lock_access_ec = threading.Lock()
lock_access_plant = threading.Lock()
lock_add_gene = threading.Lock()
lock_kegg_get = threading.Lock()

# These will hold the file path values for the programs outputs and current and project directory.
path_chem = ''
path_cwd = ''
path_fasta = ''
path_gene = ''
path_main = ''
path_raw_gene = ''

thread_lim = 5  # The max number of processor threads to be used by program.


def main():
    """
    This is the main function of the file which calls specific functions in order when running and displays the total
    run time at the end of the code execution.
    """
    init_setup()
    get_parse_pathway_genes()
    flavonoid_predictions()
    make_plant_ec_counts()
    # build_nt_fasta_by_ec()

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
    global path_cwd, path_main, path_chem, path_fasta, path_gene, list_plant_paths, list_all_plants, path_raw_gene

    # check if user supplied directory name as a commandline argument. if not, default to 'data'
    opt_dir_name = ''
    if len(sys.argv) > 1:  # alternate directory name supplied
        opt_dir_name = sys.argv[1]  # get directory name
    else:
        print('No directory name supplied, defaulting to `data`. Supply name using `python3 keggv2.py dir_name`')
        opt_dir_name = 'data'

    # update the global path values based on user decision
    path_cwd = os.getcwd() + SEP
    path_main = path_cwd + opt_dir_name
    path_chem = path_main + DIR_CHEM
    path_fasta = path_main + DIR_FASTA
    path_gene = path_main + DIR_GENE
    path_raw_gene = path_main + DIR_RAW

    # (futil.py) initialize directories if they don't exist
    init_dir(path_main)
    init_dir(path_gene)
    init_dir(path_fasta)
    init_dir(path_chem)
    init_dir(path_raw_gene)

    list_plant_paths = [i + j for i in plant_list for j in path_map_list]  # make plant and path combination list

    for key in plant_dict:  # build list of plant objects
        tmp_plant = Plant(code=key, name=plant_dict[key])  # object made for current plant
        if not tmp_plant.is_in(list_all_plants): list_all_plants.append(tmp_plant)  # add if not already present in list


def get_parse_pathway_genes():
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
        write_append(tmp_file_path, tmp_output, write_over=True)  # write the file for the pathway
    write_append(master_gene, master_output, write_over=True)  # write the master file


def path_parse(paths):
    """
    Given a pathway for a specific plant, the program then passes it to KEGG where it retrieves the appropriate entry.
    The EC number, KO number and each GENE entry are parsed from the data sent back by KEGG and then are appropriately
    saved by updating the plant objects.
    """
    for path in paths:
        global list_all_plants, list_genes_by_path
        print(path)
        no_data = False
        with lock_kegg_get:
            raw = ''
            try:
                with open(path_raw_gene + SEP + path + '.csv', 'r') as tmp_read:
                    raw = tmp_read.read()
                    tmp_read.close()
            except FileNotFoundError:
                raw = kegg.get(path)  # get KEGG entry for pathway
                with open(path_raw_gene + SEP + path + '.csv', 'w') as tmp_get:
                    try:
                        tmp_get.write(raw)
                    except TypeError:
                        no_data = True
                    tmp_get.close()

            gene_entry = kegg.parse(raw)  # parse kegg entry into dictionary for easier access
            entry_dict = gene_entry.get(GKY)  # get the data for the dictionary key GENE

        if not no_data:
            if entry_dict is not None:
                plant_code = ''.join(re.split(RE_ALPH, path))  # get only the letters from the path
                plant_name = plant_dict.get(plant_code)  # get the plant name by accessing the plant dictionary
                with lock_add_gene:
                    list_genes_by_path.append(PathGene(path=path))  # add new path to the list
                for key in entry_dict:
                    try:
                        # find EC number in the entry using regular expressions then remove square brackets
                        ec_num = re.findall(RE_EC, entry_dict[key])

                        in_count = 0
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
                                if gene_path.path == path: gene_path.genes.append(
                                    tmp_gene)  # add to list of genes for path
                        with lock_access_plant:
                            for index, plant in enumerate(list_all_plants):
                                if plant.name == tmp_gene.plant:
                                    tmp_plant = plant
                                    # add to list of plants genes if not already present
                                    if not tmp_gene.is_in(tmp_plant.genes): tmp_plant.genes.append(tmp_gene)
                                    # add to list of all genes if not already present
                                    if not tmp_gene.is_in(list_all_genes): list_all_genes.append(tmp_gene)
                                    tmp_plant.ec_nums.extend(
                                        ec_num)  # add to the plants list of ec numbers (dupes okay)
                                    list_all_plants[
                                        index] = plant  # update the list of plants with modified plant object
                    except IndexError:
                        pass  # couldn't find items using regular expression findall


def flavonoid_predictions():
    """
    This is the function that goes through each plant, looks at the list of EC numbers then applies a function in order
    to determine whether or not the plant has the required EC numbers needed to synthesize each compound.
    """

    global list_all_plants
    output_list = ''  # will have list of plants for all flavs, to be written in one file
    output_yn = ''  # will have a Y/N depending on whether or not a plant is predicted
    plant_ec_output = ''
    plant_names = 'Name\n'

    for plant in list_all_plants:
        unique_nums = []
        plant_ec_output += '\n' + plant.name + ', '
        for num in plant.ec_nums:
            if num not in unique_nums: unique_nums.append(num)
        for chem_data in data_lists:
            if predict.flav_check(getattr(predict, chem_data.code.lower()), unique_nums):
                chem_data.plants.append(plant.name)  # add plant to flavonoids list
        if unique_nums:
            for num in unique_nums: plant_ec_output += num + ', '

    for plant in list_all_plants: plant_names = plant_names + plant.name + '\n'
    write_append(path_chem + SEP + '_all-plant-names.csv', plant_names, write_over=True)

    output_yn = plant_names.replace('\n', '\t')
    # create the prediction output files for each flavonoid
    for key in data_lists:
        save_file([key.plants], key.file_name, path_chem, sep='\n')
        output_list = output_list + '\n' + key.code
        output_yn = output_yn + '\n' + key.code
        for plant in key.plants: output_list = output_list + '\t' + plant

        for item in list_all_plants:
            if item.name in key.plants:
                output_yn = output_yn + '\tY'
            else:
                output_yn = output_yn + '\tN'
        item_count = len(key.plants)
        print(key.label + ' predicted in ' + str(item_count) + ' entries.')
    write_append(path_chem + SEP + '_plant-ec-nums.csv', plant_ec_output, write_over=True)
    write_append(path_chem + SEP + '_predictions_list.csv', output_list, write_over=True)
    write_append(path_chem + SEP + '_predictions_yn.csv', output_yn, write_over=True)


def make_plant_ec_counts():
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
        if tmp_plant.has_ec_count(num):
            tmp_plant.incr_ec_count(num)  # increment count if number is already present
        else:
            # create new count object for current EC number
            tmp_count = EcCounts(number=num, count=1)
            plant.ec_counts.append(tmp_count)
    with lock_access_plant:
        list_all_plant_matrix.append(tmp_plant)  # update count matrix


def build_nt_fasta_by_ec():
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
            with request.urlopen(db_url) as db_site:
                url_data = db_site.read().decode('utf-8')
        except HTTPError or URLError as url_http_err:  # error getting html data
            print('Something went wrong with error ' + url_http_err)
            continue

        # get the header of the FASTA entry using regular expressions, &gt; is the HTML representation of >
        fasta_header = ''.join(re.findall(RE_NT_HEAD, url_data)).replace('&gt;', '>') + \
                       ' {' + plant_dict[gene.plant_code.strip()] + '}'
        fasta_body = ''.join(re.findall(RE_NT_SEQ, url_data))  # get the DNA sequence body using regular expressions
        full_fasta_entry = fasta_header + '\n' + fasta_body  # create FASTA entry string
        # create new entry object
        tmp_entry = FastaEcEntry(gene=gene.gene_id, plant=plant_dict.get(gene.plant_code), dna=full_fasta_entry)
        with lock_access_ec:
            for g in gene.ec_nums:
                tmp_ec = EcFastaCollection(ec_num=g, ec_entries=[tmp_entry])
                list_fasta_ec.append(tmp_ec)


if __name__ == '__main__':
    main()

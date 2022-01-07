# custom project library imports
from flib.fconstants import *
from flib.data_types import *
import flib.prediction_logic as predict

# other library imports
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
all_genes = []
ec_collections = []
pathgenes = []
plant_matrix = []
plant_objects = []
plant_pathways = []

# The locks are to protect the values of the global variables when using multithreading.
lock_ec = threading.Lock()
lock_gene = threading.Lock()
lock_kegg = threading.Lock()
lock_plant = threading.Lock()
lock_dbget = threading.Lock()

# These will hold the file path values for the programs outputs and current
# and project directory.
path_chem = ''
path_cwd = ''
path_fasta = ''
path_gene = ''
path_main = ''
path_raw_gene = ''
path_raw_fasta = ''

thread_lim = 5  # The max number of processor threads to be used by program.

def main():
    """
    This is the main function of the file which calls specific functions in order when
    running and displays the total run time at the end of the code execution.
    """
    setup()
    get_parse_pathway_genes()
    flavonoid_predictions()
    # make_plant_ec_counts()
    # build_nt_fasta_by_ec()
    
    runtime = datetime.datetime.now() - init_time
    print('\nRun time: ' + str(runtime))

def setup():
    """
    This is the initial setup function for the program. If the user supplies a directory
    name in the command line args, then that name will be used when outputting the data.
    If no name is supplied, then the data is outputted to the directory named data. This
    function also creates the list of pathway and plant codes based on the KEGG codes
    that can be found in the JSON files. After making the list, it creates a list of plant
    objects that will be used throughout the program.
    """
    
    # Make sure global values can be used and updated
    global path_chem
    global path_cwd
    global path_fasta
    global path_gene
    global path_main
    global path_raw_gene
    global plant_objects
    global plant_pathways
    global path_raw_fasta
    
    # Check if user supplied directory name as a commandline argument.
    # if not, default to 'data'
    opt_dir_name = ''
    if len(sys.argv) > 1:  # Alternate directory name supplied by user
        opt_dir_name = sys.argv[1]  # Get directory name
    else:
        print('No directory name supplied, defaulting to `data`. '
              'Supply name using `python3 keggv2.py dir_name`')
        opt_dir_name = 'data'
    
    # Update the global path values based on user decision
    path_cwd = os.getcwd() + SEP
    path_main = path_cwd + opt_dir_name
    path_chem = path_main + DIR_CHEM
    path_fasta = path_main + DIR_FASTA
    path_gene = path_main + DIR_GENE
    path_raw_gene = path_main + DIR_RGENE
    path_raw_fasta = path_main + DIR_RFASTA
    
    #  Initialize the output directories if they don't exist
    init_dir(path_main)
    init_dir(path_gene)
    init_dir(path_fasta)
    init_dir(path_chem)
    init_dir(path_raw_gene)
    init_dir(path_raw_fasta)
    
    # Combines plant and pathway codes.
    plant_pathways = [i + j for i in plant_list for j in path_map_list]
    
    # Go through the list of plants and create a Plant object
    # which will store related predictions.
    for key in plant_dict:
        # Call constructor to make new plant.
        tmp_plant = Plant(code=key, name=plant_dict[key])
        # Add to list if not present. Prevents duplicates.
        if not tmp_plant.is_in(plant_objects): plant_objects.append(tmp_plant)

def get_parse_pathway_genes():
    """
    This function breaks the list of plant pathways into different lists in order for
    different data to be processed at the same time using multithreading. Once all of
    the threads have finished, then the list of genes by path will be looped through
    in order to create both the gene data output files for each pathway and for the
    master file that contains all of the gene information.
    """
    global path_gene
    global plant_objects
    global plant_pathways
    global thread_lim
    
    # Chunk up the list of plant-pathway codes.
    sub_lists = list_partition(plant_pathways, thread_lim)
    threads = []  # Will be used to keep track of and kill off threads.
    
    # A thread is made for each sub-list which then passes each sub-list
    # as a parameter to the pathway parser.
    for sub_list in sub_lists:
        thread = threading.Thread(target=path_parse, args=(sub_list,))
        thread.start()
        threads.append(thread)
    
    for thread in threads: thread.join()  # Wait for each thread to die before continuing.
    
    master_output = ''  # Data for the master file which holds all gathered gene data.
    for item in pathgenes:
        tmp_file = path_gene + SEP + item.path + '.csv'
        tmp_output = ''  # Data for current gene's output file.
        for gene in item.genes:
            # Get the formatted string of gene information for the current and master files.
            tmp_output += gene.simple() + '\n'
            master_output += gene.simple() + '\n'
        write_append(tmp_file, tmp_output, write_over=True)
    write_append(path_main + SEP + 'MasterList.csv', master_output, write_over=True)

def path_parse(paths):
    """
    Given a pathway for a specific plant, the program then passes it to KEGG where it
    retrieves the appropriate entry. The EC number, KO number and each GENE entry are
    parsed from the data sent back by KEGG and then are appropriately saved by updating
    the plant objects.
    """
    for path in paths:
        global pathgenes
        global plant_objects
        
        print(path)  # Not necessary, but is nice for a gauge of progress.
        
        with lock_kegg:
            raw = ''
            
            # Look for local raw gene file to save time by preventing the call to kegg.get().
            try:
                with open(path_raw_gene + SEP + path + '.csv', 'r') as tmp_read:
                    raw = tmp_read.read()
                    tmp_read.close()
            
            # No local raw gene data file found so the data must be downloaded.
            except FileNotFoundError:
                raw = kegg.get(path)
                with open(path_raw_gene + SEP + path + '.csv', 'w') as tmp_get:
                    try: tmp_get.write(raw)
                    
                    # Error should only ever occur for plant-pathway codes that don't exist.
                    except TypeError: pass
                    tmp_get.close()
            
            gene_entry = {}
            if isinstance(raw, str):
                if len(raw) > 0:  # Don't bother parsing the empty entries.
                    kegg_entry = kegg.parse(raw)  # Parses kegg entry into dictionary.
                    gene_entry = kegg_entry.get(GKY)  # Get data from dictionary w/ key `GENE`.
        
        # If gene_entry is not empty, the statement will succeed. Otherwise ignore blanks.
        if gene_entry:
            plant_code = ''.join(re.split(RE_ALPH, path))  # Plant code is alpha-only.
            plant_name = plant_dict.get(plant_code)
            
            with lock_gene:
                pathgenes.append(PathGene(path=path))  # Make new PathGene  & add to list.
            
            for key in gene_entry:
                tmp_entry = gene_entry[key]
                try:
                    # Get all EC nums from entry using regular expressions
                    ec_nums = re.findall(RE_EC, tmp_entry)
                    
                    # Remove unwanted characters and format each EC num.
                    for i in range(0, len(ec_nums)):
                        item = multi_replace(ec_nums[i], [('[', ''), (']', ''),
                                                          (':', ''), (' ', '')])
                        ec_nums[i] = 'EC:' + item
                    
                    # Find and process the orthology ID.
                    orthology = multi_replace(quick_fetch(RE_KO, tmp_entry),
                                              [('[', ''), (']', '')])
                    
                    # Remove EC & KO in order to get compound name.
                    name = re.sub(RE_KO, '', (re.sub(RE_EC, '', tmp_entry)))
                    
                    # Call Gene constructor and pass in the parsed values.
                    tmp_gene = Gene(gene_id=key, plant=plant_name, ec_nums=ec_nums,
                                    path=path, ortho=orthology, compound=name,
                                    plant_code=plant_code)
                    with lock_gene:
                        for pathgene in pathgenes:  # Update list of genes for current path.
                            if pathgene.path == path: pathgene.genes.append(tmp_gene)
                    
                    with lock_plant:
                        for index, plant in enumerate(plant_objects):
                            if plant.name == tmp_gene.plant:
                                tmp_plant = plant
                                
                                # Add to list of plant genes if not present.
                                if not tmp_gene.is_in(tmp_plant.genes):
                                    tmp_plant.genes.append(tmp_gene)
                                
                                # Add to list of all genes if not present.
                                if not tmp_gene.is_in(all_genes):
                                    all_genes.append(tmp_gene)
                                
                                # Update the plant's EC nums. Dupes preferred
                                # (for the master count matrix).
                                tmp_plant.ec_nums.extend(ec_nums)
                                
                                # Update the plant list with the new/additional info
                                # for current plant.
                                plant_objects[index] = plant
                
                except IndexError:
                    pass  # There was nothing found using regular expressions.

def flavonoid_predictions():
    """
    This is the function that goes through each plant, looks at the list of EC numbers
    then applies a function in order to determine whether or not the plant has the
    required EC numbers needed to synthesize each compound.
    """
    
    global plant_objects
    output_list = ''  # Master list of all flavonoids and predicted plants.
    output_yn = ''  # Master list of all flavonoids & a Y/N based on predicted status
    plant_ec_output = ''  # For outputting each plant & EC numbers from their gene entries.
    plant_names = 'Name\n'  # Each plant's name, used in the first line of the Y/N output.
    
    for plant in plant_objects:
        unique_nums = []
        plant_ec_output += '\n' + plant.name + '\t'
        plant_names = plant_names + plant.name + '\n'
        
        # Go through the plant's EC nums and add each EC number to the output once.
        for num in plant.ec_nums:
            if num not in unique_nums:
                plant_ec_output += num + '\t'  # Add EC number to the output
                unique_nums.append(num)  # Ensures EC num isn't added to output again.
        
        for chem_data in flav_data_lists:  # Make the call to the prediction functions.
            if predict.flav_check(getattr(predict, chem_data.code.lower()), unique_nums):
                chem_data.plants.append(
                    plant.name)  # Add the name of the plant, if predicted.
    
    get_intersections(flav_data_lists)
    
    # Create the formatted strings for the output prediction file.
    output_yn = plant_names.replace('\n', '\t')
    for key in flav_data_lists:
        save_file([key.plants], key.file_name, path_chem, sep='\n')
        
        output_list = output_list + '\n' + key.code
        
        for plant in key.plants:
            output_list = output_list + '\t' + plant
        
        output_yn = output_yn + '\n' + key.code
        
        for item in plant_objects:
            if item.name in key.plants:
                output_yn = output_yn + '\tY'
            else:
                output_yn = output_yn + '\tN'
        
        item_count = len(key.plants)
        print(key.label + ' predicted in ' + str(item_count) + ' organisms.')
    write_append(path_chem + SEP + '_plant-ec-nums.tsv', plant_ec_output, write_over=True)
    write_append(path_chem + SEP + '_predictions_list.tsv', output_list, write_over=True)
    write_append(path_chem + SEP + '_predictions_yn.tsv', output_yn, write_over=True)

def get_intersections(fdat):
    """
    simple piece of code to find the intersections between specific flavonoids. it probably will
    be of little to no importance aside from my thesis.
    """
    agi_lu2 = set(gfo('AGI', fdat).plants).intersection(set(gfo('LU2', fdat).plants))
    write_append(path_chem + SEP + '__agi_lu2.csv', '\n'.join(agi_lu2), write_over=True)
    
    bun_hcc = set(gfo('BUN', fdat).plants).intersection(set(gfo('HCC', fdat).plants))
    write_append(path_chem + SEP + '__bun_hcc.csv', '\n'.join(bun_hcc), write_over=True)
    
    erd_nar = set(gfo('ERD', fdat).plants).intersection(set(gfo('NAR', fdat).plants))
    write_append(path_chem + SEP + '__erd_nar.csv', '\n'.join(erd_nar), write_over=True)
    
    kmp_myc_que = set(gfo('KMP', fdat).plants).intersection(set(gfo('MYC', fdat).plants),
                                                            set(gfo('QUE', fdat).plants))
    write_append(path_chem + SEP + '__kmp_myc_nar.csv', '\n'.join(kmp_myc_que), write_over=True)
    
    ec_egt_gc_kxn = set(gfo('EC', fdat).plants).intersection(set(gfo('EGT', fdat).plants),
                                                             set(gfo('GC', fdat).plants),
                                                             set(gfo('KXN', fdat).plants))
    write_append(path_chem + SEP + '__ec_egt_gc_kxn.csv', '\n'.join(ec_egt_gc_kxn), write_over=True)

def make_plant_ec_counts():
    """
    This function creates and outputs a 'matrix' relating to each species and EC number
    by running the fill_matrix function on multiple threads. For each gene entry
    containing a specific EC number, the program will increase the counter and display
    it at the end next to the appropriate EC number.
    """
    global plant_objects
    for plant in plant_objects:
        fill_count_matrix(plant)  # Update master matrix using the info from each plant.
    
    out = ''  # Output string for the master EC count matrix.
    for plant in plant_matrix:
        out += plant.name + ':\t'
        for count in plant.ec_counts:
            # Update with EC num & its occurrence.
            out += str(count.number) + ' (' + str(count.count) + ')\t'
        out += '\n'
    
    write_append(path_main + SEP + 'MasterECCountMatrix.tsv', out, write_over=True)

def fill_count_matrix(plant):
    """
    This function builds the ec counts for each list.
    """
    global plant_objects
    tmp_plant = plant
    
    for num in tmp_plant.ec_nums:
        if tmp_plant.has_ec_count(num):
            tmp_plant.incr_ec_count(num)  # Increment EC number's count.
        else:  # Need to create count object for current EC number.
            tmp_count = EcCounts(number=num, count=1)
            plant.ec_counts.append(tmp_count)
    
    with lock_plant:  # Update the count matrix.
        plant_matrix.append(tmp_plant)

def build_nt_fasta_by_ec():
    """
    This function uses multithreading and the information gathered from running
    run_path_parse in order to get the FASTA/DNA sequence for each of the gene entries
    that were found. As before, the program parses the list after all threads are done
    and then created a FASTA file for each EC number and created the Master FASTA file.
    """
    
    # Chunk up the list of all plant genes.
    sub_lists = list_partition(all_genes, thread_lim)
    threads = []  # For keeping track of each thread.
    print('getting data for ' + str(len(all_genes)) + ' genes')
    
    for sub_list in sub_lists:  # Make new thread for each sublist & pass to build_fasta().
        thread = threading.Thread(target=build_fasta, args=(sub_list,))
        thread.start()
        threads.append(thread)
    
    for thread in threads: thread.join()  # Don't continue until all threads are done.
    
    print('Starting to gather data for FASTA files...')
    master_fasta = path_main + SEP + 'MasterFASTA.fasta'
    master_output = ''
    for item in ec_collections:  # just to make sure file is cleared before populating
        item.ec_name = item.ec_name.replace('.', '-').replace(':', '')
        tmp_file_path = path_fasta + SEP + item.ec_name + '.fasta'
        write_append(tmp_file_path, '', write_over=True)
    
    for item in ec_collections:
        item.ec_name = item.ec_name.replace('.', '-').replace(':', '')
        tmp_file_path = path_fasta + SEP + item.ec_name + '.fasta'
        tmp_output = ''
        for entry in item.ec_entries:
            # add to string to be printed into specific EC file
            tmp_output += entry.simple() + '\n'
            
            # add to string to be printed into master FASTA
            master_output += entry.simple() + '\n'
        
        # write the file for current EC number
        write_append(tmp_file_path, tmp_output, write_over=False)
    
    # write the master FASTA file
    write_append(master_fasta, master_output, write_over=True)
    print('Done making the FASTA files.')

def build_fasta(genes):
    """
    This function uses the plant code and gene id in order to find the matching FASTA
    sequence using the appropriate dbget url. The pages are saved into memory as HTML
    and are parsed in order to extract the important information from the web page.
    After parsing, the FASTA sequences are added to EcFastaCollection objects in order
    to maintain proper association when writing all of the sequences out to files.
    """
    global ec_collections
    for gene in genes:
        # using the plant code and gene id create a string formatted as code:gene
        combined = gene.plant_code.strip() + ':' + \
                   gene.gene_id.replace('(RAP-DB) ', '').strip()
        
        # append code-gene string to the end of the dbget incomplete URL
        db_url = URL_DBGET + combined
        url_data = ''
        with lock_dbget:
            tmp_code = gene.plant_code
            tmp_id = gene.gene_id
            try:  # Look for local raw fasta file to save time by not having to download.
                with open(path_raw_fasta + SEP + tmp_code + tmp_id + TF, 'r') as tmp_read:
                    url_data = tmp_read.read()
                    tmp_read.close()
            
            # No local raw fasta data file found so the data must be downloaded.
            except FileNotFoundError:
                try:
                    # read the html from the dbget url
                    with request.urlopen(db_url) as db_site:
                        url_data = db_site.read().decode('utf-8')
                    
                    with open(path_raw_fasta + SEP + tmp_code + tmp_id + TF, 'w') as tmp_get:
                        tmp_get.write(url_data)
                        tmp_get.close()
                
                except (HTTPError, URLError) as url_http_err:
                    print('Something went wrong with error ' + url_http_err)
                    continue
        
        # get the header of the FASTA entry using regular expressions,
        # &gt; is the HTML representation of >
        fasta_header = ''.join(re.findall(RE_NT_HEAD, url_data)).replace('&gt;', '>') + \
                       ' {' + plant_dict[gene.plant_code.strip()] + '}'
        
        # get the DNA sequence body using regular expressions
        fasta_body = ''.join(re.findall(RE_NT_SEQ, url_data))
        
        full_fasta_entry = fasta_header + '\n' + fasta_body  # create FASTA entry string
        
        tmp_entry = FastaEcEntry(plant=plant_dict.get(gene.plant_code),
                                 gene=gene.gene_id, dna=full_fasta_entry)
        with lock_ec:
            for g in gene.ec_nums:
                tmp_ec = EcFastaCollection(ec_num=g, ec_entries=[tmp_entry])
                ec_collections.append(tmp_ec)

if __name__ == '__main__':
    main()

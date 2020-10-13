import threading
import datetime
from bioservices.kegg import KEGG
from lib.jsondata import *
from lib.miscstrings import *
from lib.datatypes import *
import sys
import re

init_time = datetime.datetime.now()
kegg = KEGG()
chem_path = ''
cwd = ''
fasta_path = ''
gene_path = ''
main_dir = ''

count_matrix = [[]]
data_lists = []
dna_list = []
enz_class_list = []
fasta_by_enz_class = []
lock_master = threading.Semaphore()
lock_parsing = threading.Semaphore()
lock_write = threading.Semaphore()
locks = []
master_gene_list = []
master_list = []
master_uniq = []
path_and_species_list = []
plant_flavs = []
species_list = full_list
thread_parsing_data = []
tmp_data_holder = []

all_plants = []
all_flavs = []

def main():
    global path_and_species_list
    path_and_species_list = [i + j for i in species_list for j in pathway_list]
    start()
    # write_readme(main_dir, fn_readme, init_time, fasta_path, gene_path)
    init_data()
    parse_main()
    apply_logic()
def start():
    global chem_path
    global fasta_path
    global gene_path
    global main_dir
    global cwd
    decision = ''
    if len(sys.argv) > 1:
        decision = sys.argv[1]
    else:
        print("No directory name supplied in args, defaulting to directory 'data'. "
              "Supply directory name in terminal using 'python3 kegg-prog.py dir_name'")
        decision = 'data'

    # create sub dirs for gene, FASTA & chemical data
    cwd = os.getcwd() + slash
    main_dir = cwd + decision
    chem_path = main_dir + chem_dir
    fasta_path = main_dir + fasta_dir
    gene_path = main_dir + gene_dir

    # replaced WindowsError with OSError for more general usage. try to make data directories and handle any errors
    try:
        os.mkdir(main_dir)
    except OSError:
        print("Error making main dir.")
        pass
    try:
        os.mkdir(gene_path)
    except OSError:
        print("Error making Gene dir.")
        pass
    try:
        os.mkdir(fasta_path)
    except OSError:
        print("Error making FASTA dir.")
        pass
    try:
        os.mkdir(chem_path)
    except OSError:
        print("Error making Chemical dir.")
        pass

def init_data():
    for li in full_list:
        np = Plant()
        np.code = li
        np.name = species_pairs.get(li)
        for pc in pathway_list:
            np.path_codes.append(li + pc)
        all_plants.append(np)
    for od in out_data:
        nf = Flavonoid()
        nf.name = od[0]
        nf.req = od[2]
        all_flavs.append(nf)

def parse_main():
    global master_list
    global master_uniq
    global locks
    print('- parsing list of species & pathways...')
    threads = []
    chunked = list(chunk(all_plants, CHUNK_SIZE))
    for plant in all_plants:
        for path in plant.path_codes:
            tmp_ecs = gene_pathway_data(path)
            for te in tmp_ecs:
                if te not in plant.ec_numbers:
                    plant.ec_numbers.append(te)
    # for plant in all_plants:
    #     print(plant.name + ' ' + plant.code + '\n' + str(plant.path_codes) + '\n' + str(plant.ec_numbers))

def gene_pathway_data(pathway_id):
    print(pathway_id)
    raw = kegg.get(pathway_id)
    gd = kegg.parse(raw)
    gene_lines = []
    fetched_genes = gd.get('GENE')
    ec_vals = []

    if fetched_genes is not None:
        gene_vals = fetched_genes.values()
        for gv in gene_vals:
            gene_lines.append(gv)
        for gene in gene_lines:
            all_ec = re.findall('(\[EC:.*\])', gene)
            ec_vals.extend(all_ec)
    return ec_vals

def apply_logic():
    for i in range(0, len(all_plants)):
        plant = all_plants[i]
        tmp_ec = plant.ec_numbers
        for j in range(0, len(all_flavs)):
            af = all_flavs[j]
            if af.req(tmp_ec):
                af.plants.append(plant.name)
                plant.flavonoids.append(af.name)
                all_flavs[j] = af
        all_plants[i] = plant

    for f in all_flavs:
        fstr = f.name + ' plants: '
        for p in f.plants:
            fstr = fstr + p + ', '
        print(fstr)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

main()

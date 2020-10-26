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

path_fasta = ' '
path_gene = ' '
path_chem = ' '
main_dir = ' '
cwd = ' '
plant_paths = []
all_genes = []
all_plants = []
all_ec_nums = []
lock_append_gene = threading.Lock()
lock_plant_rw = threading.Lock()


def main():
    init_setup()
    main_pathway_parser()


def init_setup():
    global cwd, main_dir, path_chem, path_fasta, path_gene, plant_paths, all_plants
    decision = ''
    if len(sys.argv) > 1: decision = sys.argv[1]
    else:
        print("No directory name supplied in args, defaulting to directory 'data'. "
              "Supply directory name in terminal using 'python3 kegg-prog.py dir_name'")
        decision = 'data'
    cwd = os.getcwd() + SEP
    main_dir = cwd + decision
    path_chem = main_dir + CHEM_DIR
    path_fasta = main_dir + FASTA_DIR
    path_gene = main_dir + GENE_DIR
    init_dirs(main_dir, path_gene, path_fasta, path_chem)
    plant_paths = [i + j for i in plant_list for j in path_list]
    for key in plant_dict:
        tmp_plant = Plant(code=key, name=plant_dict[key])
        all_plants.append(tmp_plant)


def main_pathway_parser():
    global plant_paths, all_plants
    for path in plant_paths:
        get_gene_data(path)
    # for plant in all_plants: print(plant.simple())


def get_gene_data(path):
    global all_genes, all_plants
    print(path)
    raw = kegg.get(path)
    gene_entry = kegg.parse(raw)
    entry_dict = gene_entry.get('GENE')
    if entry_dict is not None:
        plant_code = ''.join(re.split('[^a-zA-Z]+', path))
        plant_name = plant_dict.get(plant_code)
        # print(plant_code)
        for key in entry_dict:
            ecn = re.findall(RE_EC, entry_dict[key])
            ko = re.findall(RE_KO, entry_dict[key])
            name = re.sub(RE_EC, '', entry_dict[key])
            name = re.sub(RE_KO, '', name)

            try:
                tmp_gene = Gene(gene_id=key, plant=plant_name, ec_num=ecn[0], k_ortho=ko[0], compound=name)
                with lock_append_gene: all_genes.append(tmp_gene)
                print(tmp_gene.simple())
                with lock_plant_rw:
                    for index, plant in enumerate(all_plants):
                        if plant.name == tmp_gene.plant:
                            tmp_plant = plant
                            tmp_plant.genes.append(tmp_gene)
                            all_plants[index] = plant

            except IndexError: pass


if __name__ == '__main__':
    main()

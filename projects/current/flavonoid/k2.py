import datetime
import os
import re
import sys
import threading
import time
import warnings
import logging.config
import multiprocessing

sys.path.append(os.getcwd().replace(os.sep + 'flavonoid', ''))
import bioservices
from bioservices.kegg import KEGG
from lib.jsondata import *
from lib.datatypes import *
from lib.pathstrings import *
from lib.compoundinfo import *

init_time = datetime.datetime.now()
kegg = KEGG()
# thread_lim = multiprocessing.cpu_count()
thread_lim = 4
path_fasta = ' '
path_gene = ' '
path_chem = ' '
main_dir = ' '
cwd = ' '
plant_paths = []
all_genes = []
all_plants = []
all_ec_nums = []
genes_by_path = []
fasta_ec = []
lock_append_gene = threading.Lock()
lock_plant_rw = threading.Lock()
lock_ec = threading.Lock()


def main():
    init_setup()
    main_pathway_parser()

    prediction()

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken: ' + str(total_time))


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
    plant_chunks = list_partition(plant_paths, thread_lim)
    threads = []
    for chunk in plant_chunks:
        t = threading.Thread(target=chunk_run, args=(chunk,))
        t.start()
        threads.append(t)
    for t in threads: t.join()
    total_out = ''
    master_gene = main_dir + SEP + 'MasterList.csv'
    for gp in genes_by_path:
        tmp_file = path_gene + SEP + gp.path + CSV
        out = ''
        for gene in gp.genes:
            out += gene.simple() + '\n'
            total_out += gene.simple() + '\n'
        basic_write(tmp_file, 'w', out)
    basic_write(master_gene, 'w', total_out)

    master_fasta = main_dir + SEP + 'MasterFASTA.csv'
    fasta_out = ''
    for fec in fasta_ec:
        tmp_name = path_fasta + SEP + fec.ec_name.replace('.', '-') + CSV
        out = ''
        for entry in fec.ec_entries:
            out += entry.simple() + '\n'
            fasta_out += entry.simple() + '\n'
        basic_write(tmp_name, 'w', out)
    basic_write(master_fasta, 'w', fasta_out)

    for plant in all_plants:
        out = plant.name + ': '
        for num in plant.ec_nums: out += num + ' '
        # print(out)


def chunk_run(chunks):
    for chunk in chunks:
        get_gene_data(chunk)


def get_gene_data(path):
    time.sleep(.5)
    global all_genes, all_plants, genes_by_path
    print(path)
    raw = kegg.get(path)
    try: gene_entry = kegg.parse(raw)
    except bioservices.BioServicesError:
        print('err parsing entry ', path)
        return
    entry_dict = gene_entry.get(G_KEY)
    if entry_dict is not None:
        plant_code = ''.join(re.split(RE_ALPH, path))
        plant_name = plant_dict.get(plant_code)
        with lock_append_gene: genes_by_path.append(PathGene(path=path))
        for key in entry_dict:
            try:
                ecn = mult_replace(quick_fetch(RE_EC, entry_dict[key]), [('[', ''), (']', '')])
                ko = mult_replace(quick_fetch(RE_KO, entry_dict[key]), [('[', ''), (']', '')])
                name = re.sub(RE_KO, '', (re.sub(RE_EC, '', entry_dict[key])))
                tmp_gene = Gene(gene_id=key, plant=plant_name, ec_num=ecn, k_ortho=ko, compound=name, path=path)
                with lock_append_gene:
                    for gp in genes_by_path:
                        if gp.path == path: gp.genes.append(tmp_gene)
                with lock_plant_rw:
                    for index, plant in enumerate(all_plants):
                        if plant.name == tmp_gene.plant:
                            tmp_plant = plant
                            tmp_plant.genes.append(tmp_gene)
                            tmp_plant.ec_nums.append(ecn)
                            all_plants[index] = plant
                build_fasta(plant_code, key, ecn)
            except IndexError: pass


def build_fasta(plant_code, gene_name, ec_num):
    # time.sleep(.5)
    global fasta_ec
    # combined = plant_code + ':' + gene_name.replace('(RAP-DB) ', '')
    # raw = kegg.get(combined)
    # if raw is None: return
    # parsed = kegg.parse(raw)
    # if parsed is None: return
    # ntseq_data = parsed.get(N_KEY)
    # if ntseq_data is None: return
    # ntseq_data = ntseq_data.replace(' ', '')
    ntseq_data = 'filler string'
    tmp_entry = EntryEC(gene=gene_name, plant=plant_dict.get(plant_code), dna=ntseq_data)
    with lock_ec:
        count = 0
        for i in range(0, len(fasta_ec)):
            if fasta_ec[i].ec_name == ec_num:
                tmp = fasta_ec[i]
                tmp.ec_entries.append(tmp_entry)
                count += 1
        if count == 0:
            tmp_ec = NumEC(ec_num=ec_num, ec_entries=[tmp_entry])
            fasta_ec.append(tmp_ec)
    # if any(x.ec_name == ec_num for x in fasta_ec):
    #     pass
    # else:
    #     pass


def prediction():
    global all_plants
    for plant in all_plants:
        for chem_data in data_lists:
            if flav_check(chem_data.label, plant.ec_nums):
                chem_data.species.append(plant.name)

    for key in data_lists:
        save_file(key.species, key.file_name, path_chem)
        item_count = len(key.species)
        print(key.label + ' predicted in ' + str(item_count) + ' entries. ' +
              'Data saved in ' + path_chem + SEP + key.file_name + '.')


if __name__ == '__main__':
    main()

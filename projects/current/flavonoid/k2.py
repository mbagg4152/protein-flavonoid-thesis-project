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

init_time = datetime.datetime.now()
kegg = KEGG()
list_all_genes = []
list_all_plants = []
list_fasta_ec = []
list_gene_strings = []
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
    init_setup()
    run_path_parse()
    run_build_fasta()
    run_fill_matrix()
    prediction()

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken: ' + str(total_time))


def init_setup():
    global path_cwd, path_main, path_chem, path_fasta, path_gene, list_plant_paths, list_all_plants
    decision = ''
    if len(sys.argv) > 1: decision = sys.argv[1]
    else:
        print("No directory name supplied in args, defaulting to directory 'data'. "
              "Supply directory name in terminal using 'python3 kegg-prog.py dir_name'")
        decision = 'data'
    path_cwd = os.getcwd() + SEP
    path_main = path_cwd + decision
    path_chem = path_main + CHEM_DIR
    path_fasta = path_main + FASTA_DIR
    path_gene = path_main + GENE_DIR
    init_dirs(path_main, path_gene, path_fasta, path_chem)
    list_plant_paths = [i + j for i in plant_list for j in path_list]
    for key in plant_dict:
        tmp_plant = Plant(code=key, name=plant_dict[key])
        if not tmp_plant.is_in(list_all_plants): list_all_plants.append(tmp_plant)


def run_path_parse():
    global list_plant_paths, list_all_plants
    plant_chunks = list_partition(list_plant_paths, thread_lim)
    threads = []
    for chunk in plant_chunks:
        t = threading.Thread(target=path_parse, args=(chunk,))
        t.start()
        threads.append(t)
    for t in threads: t.join()
    total_out = ''
    master_gene = path_main + SEP + 'MasterList.csv'
    for gp in list_genes_by_path:
        tmp_file = path_gene + SEP + gp.path + CSV
        out = ''
        for gene in gp.genes:
            out += gene.simple() + '\n'
            total_out += gene.simple() + '\n'
        basic_write(tmp_file, 'w', out)
    basic_write(master_gene, 'w', total_out)

    for gene in list_all_genes:
        tmp = gene.plant_code + ':' + gene.gene_id
        if tmp not in list_gene_strings: list_gene_strings.append(tmp)


def run_build_fasta():
    chunk_genes = list_partition(list_all_genes, thread_lim)
    g_threads = []
    print('getting data for ' + str(len(list_all_genes)) + ' genes')
    for chunk in chunk_genes:
        t = threading.Thread(target=build_fasta, args=(chunk,))
        t.start()
        g_threads.append(t)

    for t in g_threads: t.join()

    master_fasta = path_main + SEP + 'MasterFASTA.csv'
    fasta_out = ''
    for fec in list_fasta_ec:
        tmp_name = path_fasta + SEP + fec.ec_name.replace('.', '-') + CSV
        out = ''
        for entry in fec.ec_entries:
            out += entry.simple() + '\n'
            fasta_out += entry.simple() + '\n'
        basic_write(tmp_name, 'w', out)
    basic_write(master_fasta, 'w', fasta_out)


def run_fill_matrix():
    global list_all_plants

    for plant in list_all_plants:
        fill_count_matrix(plant)

    out = ''
    for plant in list_all_plant_matrix:
        out += plant.name + ': '
        for count in plant.ec_counts:
            out += '[' + count.number + ' count: ' + str(count.count) + '] '
        out += '\n'
    basic_write(path_main + SEP + 'MasterCount.csv', 'w', out)


def path_parse(paths):
    # time.sleep(.5)
    for path in paths:
        global list_all_plants, list_genes_by_path
        print(path)
        with lock_kegg_get:
            raw = kegg.get(path)
            try: gene_entry = kegg.parse(raw)
            except bioservices.BioServicesError:
                print('err parsing entry ', path)
                return
            entry_dict = gene_entry.get(G_KEY)
        if entry_dict is not None:
            plant_code = ''.join(re.split(RE_ALPH, path))
            plant_name = plant_dict.get(plant_code)
            with lock_add_gene: list_genes_by_path.append(PathGene(path=path))
            for key in entry_dict:
                try:
                    ecn = mult_replace(quick_fetch(RE_EC, entry_dict[key]), [('[', ''), (']', '')])
                    if ecn not in ec_nums_of_interest: pass
                    ko = mult_replace(quick_fetch(RE_KO, entry_dict[key]), [('[', ''), (']', '')])
                    name = re.sub(RE_KO, '', (re.sub(RE_EC, '', entry_dict[key])))
                    tmp_gene = Gene(gene_id=key, plant=plant_name, ec_num=ecn, k_ortho=ko, compound=name, path=path,
                                    plant_code=plant_code)
                    with lock_add_gene:
                        for gp in list_genes_by_path:
                            if gp.path == path: gp.genes.append(tmp_gene)
                    with lock_access_plant:
                        for index, plant in enumerate(list_all_plants):
                            if plant.name == tmp_gene.plant:
                                tmp_plant = plant
                                if not tmp_gene.is_in(tmp_plant.genes): tmp_plant.genes.append(tmp_gene)
                                if not tmp_gene.is_in(list_all_genes): list_all_genes.append(tmp_gene)
                                tmp_plant.ec_nums.append(ecn)
                                list_all_plants[index] = plant
                except IndexError: pass


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

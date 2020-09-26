import threading
import datetime
from bioservices.kegg import KEGG
from lib.json_data import *
from lib.misc_strings import *
from lib.project_classes import *

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


def main():
    global path_and_species_list
    path_and_species_list = [i + j for i in species_list for j in pathway_list]
    start()
    parse_master()
    make_matrix_and_counts()
    make_fasta()
    write_readme(main_dir, fn_readme, init_time, fasta_path, gene_path)
    finish_up()


def start():
    global chem_path
    global fasta_path
    global gene_path
    global main_dir
    global cwd
    # check if main dir exists
    # if os.path.exists(main_dir):
    #     decision = input("Directory already exists. Press return key to continue, or type new name: ")
    #     # if return key hit will continue running th code, will overwrite anything in the preexisting folder
    #     if decision == "":
    #         pass
    #     else:
    #         # stops code completely and the code will need to be restarted with a different name
    #         sys.exit("Try an unused folder name next time")

    # create sub dirs for gene, FASTA & chemical data
    cwd = os.getcwd() + slash
    main_dir = cwd + 'data'
    chem_path = main_dir + chem_dir
    fasta_path = main_dir + fasta_dir
    gene_path = main_dir + gene_dir

    # replaced WindowsError with OSError for more general usage
    # try to make data directories and handle any errors
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


# function that fetches the required data
def gene_pathway_data(pathway_id):
    print(pathway_id)
    entry_lines = str(kegg.get(pathway_id)).split(NL)  # gets all of the data and splits it by line
    # print genes
    line_count = 0
    gene_locator = 0
    compound_locator = 0
    for line in entry_lines:  # finds the places that have the genes listed
        entry_lines[line_count] = line.strip()  # remove unwanted whitespace
        if line.startswith('GENE'):  # finds where GENE is at in each entry
            entry_lines[line_count] = line.replace("GENE", "").strip()
            gene_locator = line_count  # gene locator is now the item in the list that has GENE
        if line.startswith('COMPOUND'):  # finds where COMPOUND is in the entry
            compound_locator = line_count  # now the item in the list that begins with compound
        line_count += 1
    gene_lines = entry_lines[gene_locator:compound_locator]  # makes a list that is just the gene entry lines
    line_count = 0
    for gene in gene_lines:  # this section makes a list of lists that are appreciatively separated
        # makes ^*^ the signifier for splitting
        split_sig = '^*^'
        gene = gene.replace('  ', split_sig).replace(';', split_sig).replace('[', split_sig).replace(']', '')
        gene = gene.split(split_sig)
        alpha_only = NIX
        for char in pathway_id:
            if char.isalpha():
                alpha_only += char

        gene.insert(0, species_pairs[alpha_only])
        j_count = 0
        for g in gene:  # cleans up list of lists of extra spaces at the beginning and end of each list
            gene[j_count] = g.strip()
            j_count += 1
        gene_lines[line_count] = gene  # replaces the list with the new cleaned list of lists
        # print(str(gene))
        line_count += 1  # iterates through each list in the entry
    return gene_lines


def parse_helper(plant_paths, sem_index):
    global master_list
    global tmp_data_holder
    global thread_parsing_data

    for path in plant_paths:
        lock_parsing.acquire()
        try:
            tmp_data_holder[sem_index] = gene_pathway_data(path)
        finally:
            lock_parsing.release()
        lock_write.acquire()
        try:
            thread_parsing_data[sem_index].extend(tmp_data_holder[sem_index])
            try:
                save_file(tmp_data_holder[sem_index], GDATA + path + CSV, gene_path)
            except AttributeError:
                pass
        finally:
            lock_write.release()

    lock_master.acquire()
    try:
        master_list.extend(thread_parsing_data[sem_index])
    finally:
        lock_master.release()


def parse_master():
    global master_list
    global master_uniq
    global locks
    print('- parsing list of species & pathways...')
    threads = []
    chunked = list(chunk(path_and_species_list, CHUNK_SIZE))
    for i in range(0, len(chunked)):
        locks.append(threading.Semaphore(1))
        thread_parsing_data.append([])
        tmp_data_holder.append([])
        # parse_sems.append(threading.Semaphore())

    chunk_index = 0
    for chunks in chunked:
        t = threading.Thread(target=parse_helper, args=(chunks, chunk_index))
        t.start()
        threads.append(t)
        chunk_index += 1

    for t in threads:
        t.join()
    master_uniq = remove_dupes(master_list)
    master_uniq = list(filter(None, master_uniq))
    count = 0
    for i in master_uniq:  # removes false values and iterates through the list of lists
        master_uniq[count] = list(filter(None, master_uniq[count]))
        count += 1


def make_matrix_and_counts():
    global master_gene_list
    global count_matrix

    ec_list = unique_element_list(master_uniq, 'last')

    # creating the matrix and adding up the counts
    count_matrix = [['Species']]
    # adds the Unique EC numbers to the end of the matrix
    count_matrix[0].extend(ec_list)
    for i in species_list:  # @#
        # first item in each row (but first) is the matrix
        count_matrix.append([i])

    cols = 0
    for outer in count_matrix[0]:
        if cols != 0:  # first column isn't actually an EC#,
            ec = outer
            rows = 0
            for inner in count_matrix:  # for each species (using the species code)
                if rows != 0:  # first row isn't actually a species
                    species = species_pairs[inner[0]]
                    counter3 = 0
                    for unique in master_uniq:  # iterate over the culled master list to check for matching sets
                        if unique[0] == species and unique[len(unique) - 1] == ec:
                            counter3 += 1
                    count_matrix[rows].append(str(counter3))
                rows += 1
        cols += 1

    # change master count to be actual species:
    count = 0
    for i in count_matrix:
        if count != 0:
            count_matrix[count][0] = species_pairs[i[0]]  # replaces species code with genus specie names
        count += 1

    # Make Master Files
    print('- making master files: no dupes & count matrix')
    save_file(master_uniq, 'Master_List.csv', main_dir)
    save_file(count_matrix, 'Master_Count.csv', main_dir)

    # make a master fasta file
    print('- about to make master FASTA')
    swapped_order = {v: k for k, v in species_pairs.items()}  # reverses dictionary keys and values
    master_gene_list = []
    for i in master_uniq:
        # combines species codes and gene numbers in a list to be used for the master fasta function
        master_gene_list.append(swapped_order[i[0]] + ':' + i[1])


def get_master_fasta(gene):
    print('- fetching data for master FASTA...')
    global dna_list, master_gene_list
    threads = []
    chunked = list(chunk(master_gene_list, CHUNK_SIZE))
    for chunks in chunked:
        try:
            thread = threading.Thread(target=master_helper, args=(chunks,))
            thread.start()
            threads.append(thread)

        except:
            pass

    for thread in threads:
        thread.join()
    return dna_list


def master_helper(gene_chunk):
    for gene in gene_chunk:
        # calls the entry from KEGG and splits it into new lines
        gene_fasta_data = str(kegg.get(gene)).split(NL)
        global dna_list
        line_count = 0
        ntseq_locator = 0
        organism_name = ''
        for line in gene_fasta_data:
            # removes blank spaces at the beginning and end of each line
            gene_fasta_data[line_count] = line.strip()
            if line.startswith('ORGANISM'):  # finds where the entry that begins with organism is
                gene_fasta_data[line_count] = line.replace("ORGANISM", "").strip()
                find_blanks = gene_fasta_data[line_count].find(" ")
                organism_name = gene_fasta_data[line_count][find_blanks:]
            line_count += 1
        line_count = 0
        for line in gene_fasta_data:
            gene_fasta_data[line_count] = line.strip()
            if line.startswith('ORTHOLOGY'):
                gene_fasta_data[line_count] = line.strip
                find_ec = line.find("EC:")
                gene_fasta_data[line_count] = line[find_ec:-1].replace("[", "").replace("EC:", "")
                ec_number = "EC " + gene_fasta_data[line_count]  # adds EC back
                # adds > to beginning to find the beginning of each entry more easily, adds % between EC number and
                # the rest of the entry to help separate the EC number for later and removes semicolon from the gene
                # entry and adds the gene number
                joined_organism_ec = [">" + str(organism_name).strip() + "%" + ec_number + "%" + gene.split(":")[1]]
                dna_list.append(joined_organism_ec)  # adds the entry to the blank list
            line_count += 1
        line_count = 0
        for line in gene_fasta_data:
            gene_fasta_data[line_count] = line.strip()
            if line.startswith('NTSEQ'):
                gene_fasta_data[line_count] = line.replace("NTSEQ", "").strip()
                ntseq_locator = line_count
            line_count += 1
        dna_data_list = gene_fasta_data[ntseq_locator:]
        dna_seq = dna_data_list[1:len(dna_data_list) - 2]  # Takes just the DNA sequence
        sep = ''
        # combines the separate DNA sequence lines into one string and turns that into a single entry list
        joined_dna_seq = [sep.join(dna_seq)]
        dna_list.append(joined_dna_seq)  # adds single entry list to the list of lists


def make_fasta():
    print('- saving master fasta...')
    master_fasta = get_master_fasta(master_gene_list)
    save_file(master_fasta, 'Master_FASTA.csv', fasta_path)

    # master_fasta = get_master_fasta(gene_list_from_master)  # should actually go above the first time it gets called
    # make a fasta file for each EC number saved in the FASTA dir

    counter = 0
    print('- looping through master fasta...')
    for i in master_fasta:
        if i[0].startswith('>'):
            split_ec = i[0].split('%')  # finds the EC numbers using the % added previously

            this_ec = split_ec[1]
            ec_count = 0
            ec_flag = 'false'
            counter2 = 0
            for j in enz_class_list:  # if the EC number is already present sets it to true and continues
                if this_ec == j:
                    ec_count = counter2
                    ec_flag = 'true'
                counter2 += 1
            if ec_flag == 'false':  # if the EC number is not there, it will be added to the list
                enz_class_list.append(this_ec)
                fasta_by_enz_class.append(master_fasta[counter:counter + 2])
            else:
                fasta_by_enz_class[ec_count].extend(master_fasta[counter:counter + 2])

        counter += 1

    # creates the FASTA files by EC numbers
    print('- creating fasta files by ec number...')
    counter = 0
    for i in enz_class_list:
        name = i.replace('.', 'p').replace('EC ', NIX).replace(SP, NIX)
        save_file(fasta_by_enz_class[counter], name + CSV, fasta_path)
        # print('- Found data for ' + name)
        counter += 1


def finish_up():
    global data_lists

    master_ec_list = [['species', 'EC#s']]
    counter = 0
    print('- filling master matrix...')
    for i in count_matrix:  # for each species
        species_ec_list = []
        if counter == 0:
            pass

        else:
            counter2 = 0
            for j in i:  # for EC in species
                if counter2 == 0:
                    species_ec_list.append(j)
                else:
                    if str(j) == "0":
                        pass
                    else:
                        species_ec_list.append(count_matrix[0][counter2])
                counter2 += 1
        master_ec_list.append(species_ec_list)
        counter += 1

    master_ec_list = master_ec_list[1:]

    for item in out_data:
        tmp_data = ChemData(item[0], item[1], item[2], item[3])
        data_lists.append(tmp_data)

    print('- looping through master ec list...')
    for i in master_ec_list:
        print(str(i))
        if len(i) > 0:
            tmp_entry = Species(i[0], 0, [])
            for chem_data in data_lists:
                if chem_data.logic in i:
                    chem_data.species.append([i[0]])
                    tmp_entry.flavonoids.append(chem_data.label)
                    tmp_entry.count += 1
            plant_flavs.append(tmp_entry)

    for key in data_lists:
        save_file(key.species, key.file_name, chem_path)
        print('\n' + key.label + ':')
        out = ''
        for li in key.species:
            out = out + str(*li) + ' | '
        print(out)

    for pf in plant_flavs:
        data_out = pf.species_string()
        print(data_out)

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken: ' + str(total_time))


main()

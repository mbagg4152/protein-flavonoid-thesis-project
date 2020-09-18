import sys
import threading

from bioservices.kegg import KEGG
import datetime
import time
from defs import *

init_time = datetime.datetime.now()
# changed from k -> kegg. 1 letter var names should only be used in iteration
kegg = KEGG()

# get active directory
cwd = os.getcwd() + slash

# get dir name from user & append to current directory
# main_dir = cwd + input("Input a save folder name: ")
main_dir = cwd + 'data'
chem_path = ''
fasta_path = ''
gene_path = ''


def start():
    global chem_path
    global fasta_path
    global gene_path
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


start()


# removes duplicate elements
def remove_dupes(dupe_list):
    unique_list = []  # creates an empty list
    for item in dupe_list:
        if item not in unique_list:
            unique_list.append(item)  # adds item to empty list if it's not already in the list
    return unique_list


# change which list should be used for test input
species_list = full_list

# this is the full list of every pathway and species from both lists
path_and_species_list = [i + j for i in species_list for j in pathway_list]

gene_list_from_master = []
ec_order_list = []
fasta_by_ec = []
master_list = []
master_no_dupes = []
master_count_matrix = [[]]


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
        line_count += 1  # iterates through each list in the entry
    return gene_lines


# function that saves each file with a name that includes pathwayID  using the data from genes_lines
# list_to_write should be a list of lists, output_dir should include an appropriate extension
def save_file(lists_to_write, output_dir, current):
    os.chdir(current)
    form_counter = 0
    # open file to be written
    writedoc = open(output_dir, 'w')
    for line in lists_to_write:
        for item in line:
            # removes the new lines in each list of list
            item = str(item).replace(NL, NIX)
            if item == NIX:  # if the list in the list of list is empty writes a dash
                writedoc.write('-')
            else:  # write the entry in the list of lists to the file
                writedoc.write(item)
            writedoc.write(', ')  # tab delineated; use "," for csv files
            if form_counter == 7:
                form_counter = 0
                writedoc.write(NL)
            form_counter += 1
    writedoc.write(NL)
    writedoc.close()


def parse_master():
    global master_list
    global master_no_dupes
    print('- parsing list of species & pathways...')
    threads = []

    for path_id in path_and_species_list:
        not_present = 0
        current_list = []
        t = threading.Thread(target=gene_pathway_data, args=(path_id,))
        try:  # need to ignore everything if there is no pathway for that species
            t.start()
            threads.append(t)
            current_list = gene_pathway_data(path_id)
        except AttributeError:
            # print("- No data found in " + path_id)
            not_present = 1
        if not_present == 0:
            master_list.extend(current_list)

        # run the function that saves each file
        try:  # try actually does it if it works
            # save_file(gene_pathway_data(path_id), GDATA + path_id + CSV, gene_path)
            save_file(current_list, GDATA + path_id + CSV, gene_path)
        except AttributeError:
            pass
    for t in threads:
        t.join()
    master_no_dupes = remove_dupes(master_list)
    # removes false values and turns it into a list
    master_no_dupes = list(filter(None, master_no_dupes))
    count = 0
    for i in master_no_dupes:  # removes false values and iterates through the list of lists
        master_no_dupes[count] = list(filter(None, master_no_dupes[count]))
        count += 1


parse_master()


# find unique EC numbers so have a generic function and run it
# be careful as there are cases of one less item - use "last" to fix that problem here
def unique_element_list(list_name, index):
    original_index = index
    element_list = []
    for i in list_name:
        # print
        if original_index == 'last':
            # assigns the string "last" to the very last list in the list of lists
            index = len(i) - 1
        if i[int(index)] not in element_list:  # finds unique EC number not in the list
            element_list.append(i[int(index)])  # adds it to the list
    return element_list


def make_matrix_and_counts():
    global gene_list_from_master
    global master_count_matrix

    ec_list = unique_element_list(master_no_dupes, 'last')

    # creating the matrix and adding up the counts
    master_count_matrix = [['Species']]
    # adds the Unique EC numbers to the end of the matrix
    master_count_matrix[0].extend(ec_list)
    for i in species_list:  # @#
        # first item in each row (but first) is the matrix
        master_count_matrix.append([i])

    cols = 0
    for outer in master_count_matrix[0]:
        if cols != 0:  # first column isn't actually an EC#,
            ec = outer
            rows = 0
            for inner in master_count_matrix:  # for each species (using the species code)
                if rows != 0:  # first row isn't actually a species
                    species = species_pairs[inner[0]]
                    counter3 = 0
                    for unique in master_no_dupes:  # iterate over the culled master list to check for matching sets
                        if unique[0] == species and unique[len(unique) - 1] == ec:
                            counter3 += 1
                    master_count_matrix[rows].append(str(counter3))
                rows += 1
        cols += 1

    # change master count to be actual species:
    count = 0
    for i in master_count_matrix:
        if count != 0:
            master_count_matrix[count][0] = species_pairs[i[0]]  # replaces species code with genus specie names
        count += 1

    # Make Master Files
    print('- making master files: no dupes & count matrix')
    save_file(master_no_dupes, 'Master_List.csv', main_dir)
    save_file(master_count_matrix, 'Master_Count.csv', main_dir)

    # make a master fasta file
    print('- about to make master FASTA')
    swapped_order = {v: k for k, v in species_pairs.items()}  # reverses dictionary keys and values
    gene_list_from_master = []
    for i in master_no_dupes:
        # combines species codes and gene numbers in a list to be used for the master fasta function
        gene_list_from_master.append(swapped_order[i[0]] + ':' + i[1])


make_matrix_and_counts()

dna_info_list = []


def chunk(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_master_fasta(gene):
    print('- fetching data for master FASTA...')
    global dna_info_list, gene_list_from_master
    threads = []
    chunked = list(chunk(gene_list_from_master, CHUNK_SIZE))
    for chunks in chunked:
        try:
            thread = threading.Thread(target=master_helper, args=(chunks,))
            thread.start()
            threads.append(thread)

        except:
            pass

    for thread in threads:
        thread.join()
    return dna_info_list


def master_helper(gene_chunk):
    for gene in gene_chunk:
        # calls the entry from KEGG and splits it into new lines
        gene_fasta_data = str(kegg.get(gene)).split(NL)
        global dna_info_list
        lines = 0
        ntseq_locator = 0
        organism_name = ''
        for line in gene_fasta_data:
            # removes blank spaces at the beginning and end of each line
            gene_fasta_data[lines] = line.strip()
            if line.startswith('ORGANISM'):  # finds where the entry that begins with organism is
                gene_fasta_data[lines] = line.replace("ORGANISM", "").strip()
                find_blanks = gene_fasta_data[lines].find(" ")
                organism_name = gene_fasta_data[lines][find_blanks:]
            lines += 1
        lines = 0
        for line in gene_fasta_data:
            gene_fasta_data[lines] = line.strip()
            if line.startswith('ORTHOLOGY'):
                gene_fasta_data[lines] = line.strip
                find_ec = line.find("EC:")
                gene_fasta_data[lines] = line[find_ec:-1].replace("[", "").replace("EC:", "")
                ec_number = "EC " + gene_fasta_data[lines]  # adds EC back
                # adds > to beginning to find the beginning of each entry more easily, adds % between EC number and
                # the rest of the entry to help separate the EC number for later and removes semicolon from the gene
                # entry and adds the gene number
                joined_organism_ec = [">" + str(organism_name).strip() + "%" + ec_number + "%" + gene.split(":")[1]]
                dna_info_list.append(joined_organism_ec)  # adds the entry to the blank list
            lines += 1
        lines = 0
        for line in gene_fasta_data:
            gene_fasta_data[lines] = line.strip()
            if line.startswith('NTSEQ'):
                gene_fasta_data[lines] = line.replace("NTSEQ", "").strip()
                ntseq_locator = lines
            lines += 1
        dna_data_list = gene_fasta_data[ntseq_locator:]
        dna_seq = dna_data_list[1:len(dna_data_list) - 2]  # Takes just the DNA sequence
        sep = ''
        # combines the separate DNA sequence lines into one string and turns that into a single entry list
        joined_dna_seq = [sep.join(dna_seq)]
        dna_info_list.append(joined_dna_seq)  # adds single entry list to the list of lists


def make_fasta():
    print('- saving master fasta...')
    master_fasta = get_master_fasta(gene_list_from_master)
    save_file(get_master_fasta(gene_list_from_master), 'Master_FASTA.csv', fasta_path)

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
            for j in ec_order_list:  # if the EC number is already present sets it to true and continues
                if this_ec == j:
                    ec_count = counter2
                    ec_flag = 'true'
                counter2 += 1
            if ec_flag == 'false':  # if the EC number is not there, it will be added to the list
                ec_order_list.append(this_ec)
                fasta_by_ec.append(master_fasta[counter:counter + 2])
            else:
                fasta_by_ec[ec_count].extend(master_fasta[counter:counter + 2])

        counter += 1

    # creates the FASTA files by EC numbers
    print('- creating fasta files by ec number...')
    counter = 0
    for i in ec_order_list:
        name = i.replace('.', 'p').replace('EC ', NIX).replace(SP, NIX)
        save_file(fasta_by_ec[counter], name + CSV, fasta_path)
        # print('- Found data for ' + name)
        counter += 1


make_fasta()


def write_readme():
    # Creates ReadMe file
    print('- creating README...')
    with open(main_dir + readme, 'w') as readme_doc:
        readme_doc.write("KEGG_v1p1.py\n")
        readme_doc.write(init_time.strftime("%m-%d-%Y") + "\n")
        readme_doc.write(main_dir + "\n")
        readme_doc.write(
            'This script creates a series of files related to the genes associated with plant flavonoids ' +
            'from various species of plants. This script first creates the MasterCount and MasterList ' +
            'files; the MasterCount counts the number genes each plant species have that correspond with ' +
            'each EC number; while the MasterList lists every gene with number for each plant specie. ' +
            'These are located in ' + os.getcwd() + '. The script also creates files that only contains ' +
            'the genes of a single plant species biochemical pathway which are located in ' + gene_path +
            '. The script also creates a Master FASTA files which contains the DNA sequence of each gene' +
            ' and FASTA files organized by EC number, these are located in ' + fasta_path)
        readme_doc.close()


write_readme()

masterEC_list = [['species', 'EC#s']]
counter = 0
print('- filling master matrix...')
for i in master_count_matrix:  # for each species
    species_epic_list = []
    if counter == 0:
        pass

    else:
        counter2 = 0
        for j in i:  # for EC in species
            if counter2 == 0:
                species_epic_list.append(j)
            else:
                if str(j) == "0":
                    pass
                else:
                    species_epic_list.append(master_count_matrix[0][counter2])
            counter2 += 1
    masterEC_list.append(species_epic_list)
    counter += 1

masterEC_list = masterEC_list[1:]

out_data = [(APIG, [], apigenin, apig_file), (BUTE, [], butein, bute_file), (CATE, [], catechin, cate_file),
            (CYAN, [], cyanidin, cyan_file), (ECAT, [], epicatechin, ecat_file),
            (EPIG, [], epigallocatechin, epig_file), (ERIO, [], eriodictyol, erio_file),
            (GALL, [], gallocatechin, gall_file), (GENI, [], genistein, geni_file),
            (KAEM, [], kaempferol, kaem_file), (LUTE, [], luteolin, lute_file), (MYRI, [], myricetin, myri_file),
            (NARI, [], naringenin, nari_file), (QUER, [], quercetin, quer_file), (EC1, [], EC1, ec1_file),
            (EC2, [], EC2, ec2_file), (EC3, [], EC3, ec3_file), (EC4, [], EC4, ec4_file)]

data_lists = []
spec_flavs = []
for item in out_data:
    tmp_data = ChemData(item[0], item[1], item[2], item[3])
    data_lists.append(tmp_data)

print('- looping through master ec list...')
for i in masterEC_list:
    if len(i) > 0:
        tmp_entry = Species(i[0], 0, [])
        for c_data in data_lists:
            if c_data.logic in i:
                c_data.species.append([i[0]])
                print('species = ' + str(i[0]) + 'data label = ' + c_data.label)
                tmp_entry.flavonoids.append(c_data.label)
                tmp_entry.count += 1
        spec_flavs.append(tmp_entry)

for key in data_lists:
    save_file(key.species, key.file_name, chem_path)
    print('\n' + key.label + ':')
    out = ''
    for li in key.species:
        out = out + str(*li) + ', '
    print(out)

    # spec_output = ''
    # for spec in spec_flavs:
    #     print(spec.species_string())
    #     spec_output += spec.species_string() + '\n'

# os.chdir(main_dir + chem_dir)
# new_file = open(slash + 'spec_file.txt', 'x')
# new_file.close()
# spec_file = open(slash + 'spec_file.txt', 'w')
# spec_file.write(spec_output)
# spec_file.close()

end_time = datetime.datetime.now()
total_time = end_time - init_time
print('\ntotal time taken: ' + str(total_time))

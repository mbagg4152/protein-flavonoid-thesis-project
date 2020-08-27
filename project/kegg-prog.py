import sys

from bioservices.kegg import KEGG

from defs import *

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
    print('- working on appending gene data...')
    entry_lines = kegg.get(pathway_id).split(NL)  # gets all of the data and splits it by line
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
            writedoc.write(',')  # tab delineated; use "," for csv files
        writedoc.write(NL)
    writedoc.write(NL)
    writedoc.close()
    # print('- Saved ' + output_dir)


def parse_master():
    global master_list
    global master_no_dupes
    print('- parsing list of species & pathways...')
    for path_id in path_and_species_list:
        not_present = 0
        current_list = []
        try:  # need to ignore everything if there is no pathway for that species
            current_list = gene_pathway_data(path_id)
        except AttributeError:
            # print("- No data found in " + path_id)
            not_present = 1
        if not_present == 0:
            master_list.extend(current_list)

        # run the function that saves each file
        try:  # try actually does it if it works
            save_file(gene_pathway_data(path_id), GDATA + path_id + CSV, gene_path)
        except AttributeError:
            pass

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
    global gene_list_from_master, master_count_matrix

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
    rev_dict = {v: k for k, v in species_pairs.items()}  # reverses dictionary keys and values
    gene_list_from_master = []
    for i in master_no_dupes:
        # combines species codes and gene numbers in a list to be used for the master fasta function
        gene_list_from_master.append(rev_dict[i[0]] + ':' + i[1])


make_matrix_and_counts()


def get_master_fasta(gene):
    print('- fetching data for master FASTA...')
    dna_info_list = []
    global gene_list_from_master
    for gene in gene_list_from_master:
        # calls the entry from KEGG and splits it into new lines
        gene_fasta_data = kegg.get(gene).split(NL)
        # print('got gene fasta data')
        lines = 0
        ntseq_locator = 0
        organism_name = ''
        # print('looking for ORGANISM in each line of gene fasta data')
        for line in gene_fasta_data:
            # removes blank spaces at the beginning and end of each line
            gene_fasta_data[lines] = line.strip()
            if line.startswith('ORGANISM'):  # finds where the entry that begins with organism is
                # removes organism and removes blank spaces ta beginning and end
                gene_fasta_data[lines] = line.replace("ORGANISM", "").strip()
                find_blanks = gene_fasta_data[lines].find(
                    " ")  # finds where the double blank space is in the organism line
                organism_name = gene_fasta_data[lines][
                                find_blanks:]  # the genus specie name is left from the orginal entry
            lines += 1
        lines = 0
        # print('looking for ORTHOLOGY in each line of gene fasta data')
        for line in gene_fasta_data:
            gene_fasta_data[lines] = line.strip()
            if line.startswith('ORTHOLOGY'):
                gene_fasta_data[lines] = line.strip
                # finds within the line where the EC number is
                find_ec = line.find("EC:")
                # removes the beginning bracket and EC: leaving just the number
                gene_fasta_data[lines] = line[find_ec:-1].replace("[", "").replace("EC:", "")
                ec_number = "EC " + gene_fasta_data[lines]  # adds EC back
                # adds > to beginning to find the beginning of each entry more easily, adds % between EC number and
                # the rest of the entry to help separate the EC number for later and removes semicolon from the gene
                # entry and adds the gene number
                joined_organism_ec = [">" + str(organism_name).strip() + "%" + ec_number + "%" + gene.split(":")[1]]
                dna_info_list.append(joined_organism_ec)  # adds the entry to the blank list
            lines += 1
        lines = 0
        # print('looking for NTSEQ in each line of gene fasta data')
        for line in gene_fasta_data:
            gene_fasta_data[lines] = line.strip()
            if line.startswith('NTSEQ'):
                gene_fasta_data[lines] = line.replace("NTSEQ", "").strip()
                ntseq_locator = lines
            lines += 1
        dna_data_list = gene_fasta_data[ntseq_locator:]
        # print('- len of dna data list: ' + str(len(dna_data_list)))
        dna_seq = dna_data_list[1:len(dna_data_list) - 2]  # Takes just the DNA sequence
        sep = ''
        # combines the separate DNA sequence lines into one string and turns that into a single entry list
        joined_dna_seq = [sep.join(dna_seq)]
        # print('- joined dna: ' + str(joined_dna_seq))
        # print('- gene list from master len: ' + str(len(gene_list_from_master)))
        # print('- made list of joined dna sequences. joined len: ' + str(len(joined_dna_seq)))
        dna_info_list.append(joined_dna_seq)  # adds single entry list to the list of lists
        # print('- added joined sequence list to dna info. dna info len: ' + str(len(dna_info_list)))
    return dna_info_list


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

            EC = split_ec[1]
            ec_count = 0
            ec_flag = 'false'
            counter2 = 0
            for j in ec_order_list:  # if the EC number is already present sets it to true and continues
                if EC == j:
                    ec_count = counter2
                    ec_flag = 'true'
                counter2 += 1
            if ec_flag == 'false':  # if the EC number is not there, it will be added to the list
                ec_order_list.append(EC)
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
        readme_doc.write(NOW.strftime("%m-%d-%Y") + "\n")
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
    # print "species :", i[0]
    species_epic_list = []
    if counter == 0:
        pass

    else:
        counter2 = 0
        for j in i:  # for EC in species
            if counter2 == 0:
                species_epic_list.append(j)
                # print(j)
            else:
                if str(j) == "0":
                    # print(" " + str(j) + " = no enzyme")
                    pass
                else:
                    # print(master_count_matrix[0][counter2] + " = " + str(j))
                    species_epic_list.append(master_count_matrix[0][counter2])
            counter2 += 1
    masterEC_list.append(species_epic_list)
    counter += 1

# print(masterEC_list)
masterEC_list = masterEC_list[1:]

apig_list = []
bute_list = []
cate_list = []
cyan_list = []
ecat_list = []
epig_list = []
erio_list = []
gall_list = []
kaem_list = []
lute_list = []
nari_list = []
quer_list = []
# Be careful in making these of parentheses
print('- looping through master ec list...')
for i in masterEC_list:
    if epicatechin in i:
        ecat_list.append([i[0]])
    if catechin in i:
        cate_list.append([i[0]])
    if eriodictyol in i:
        erio_list.append([i[0]])
    if apigenin in i:
        apig_list.append([i[0]])
    if luteolin in i:
        lute_list.append([i[0]])
    if naringenin in i:
        nari_list.append([i[0]])
    if butein in i:
        bute_list.append([i[0]])

    if kaempferol in i:
        kaem_list.append([i[0]])
    if quercetin in i:
        quer_list.append([i[0]])
    if cyanidin in i:
        cyan_list.append([i[0]])
    if epigallocatechin in i:
        epig_list.append([i[0]])
    if gallocatechin in i:
        gall_list.append([i[0]])
save_file(ecat_list, ecat_file, chem_path)
save_file(cate_list, cate_file, chem_path)
save_file(erio_list, erio_file, chem_path)
save_file(lute_list, lute_file, chem_path)
save_file(nari_list, nari_file, chem_path)
save_file(bute_list, bute_file, chem_path)
save_file(apig_list, apig_file, chem_path)
save_file(kaem_list, kaem_file, chem_path)
save_file(quer_list, quer_file, chem_path)
save_file(cyan_list, cyan_file, chem_path)
save_file(epig_list, epig_file, chem_path)
save_file(gall_list, gall_file, chem_path)
print('ecat data:\n\t' + str(ecat_list))
print('cate data:\n\t' + str(cate_list))
print('erio data:\n\t' + str(erio_list))
print('nari data:\n\t' + str(nari_list))
print('lute data:\n\t' + str(lute_list))
print('bute data:\n\t' + str(bute_list))
print('apig data:\n\t' + str(apig_list))
print('kaem data:\n\t' + str(kaem_list))
print('quer data:\n\t' + str(quer_list))
print('cyan data:\n\t' + str(cyan_list))
print('epig data:\n\t' + str(epig_list))
print('gall data:\n\t' + str(gall_list))

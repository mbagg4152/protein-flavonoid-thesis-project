import json
import os
from lib.miscvals import *
from lib.pathstrings import SEP

def get_json_data(file_name, key):
    data = ''
    try:
        with open(file_name) as jsonFile:
            data = json.load(jsonFile)
    except FileNotFoundError:
        file_name = '..' + SEP + file_name
        with open(file_name) as jsonFile:
            data = json.load(jsonFile)
    return data[key]
# removes duplicate elements
def remove_dupes(dupe_list):
    unique_list = []  # creates an empty list
    for item in dupe_list:
        # adds item to empty list if it's not already in the list
        if item not in unique_list: unique_list.append(item)
    return unique_list

# find unique EC numbers so have a generic function and run it
# be careful as there are cases of one less item - use "last" to fix that problem here
def unique_element_list(list_name, index):
    original_index = index
    element_list = []
    for i in list_name:  # assigns the string "last" to the very last list in the list of lists
        if original_index == 'last': index = len(i) - 1  # finds unique EC num not in the list & adds it to the list
        if i[int(index)] not in element_list: element_list.append(i[int(index)])
    return element_list

def list_partition(items, number):
    # looping till length l
    for i in range(0, len(items), number): yield items[i:i + number]

def write_readme(main_dir, readme, init_time, fasta_path, gene_path):
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

def save_file(lists_to_write, output_dir, current):
    os.chdir(current)
    form_counter = 0
    # open file to be written
    writedoc = open(output_dir, 'w')
    for line in lists_to_write:
        for item in line:
            item = str(item).replace(NL, NIX)  # removes the new lines in each list of list
            if item == NIX:   writedoc.write('-')  # if the list in the list of list is empty writes a dash
            else: writedoc.write(item)  # write the entry in the list of lists to the file
            writedoc.write(', ')  # tab delineated; use "," for csv files
            if form_counter == 7:
                form_counter = 0
                writedoc.write(NL)
            form_counter += 1
    writedoc.write(NL)
    writedoc.close()

def basic_write(path, mode, content):
    file = open(path, mode)
    file.write(content)
    file.close()

def is_http_error(msg):
    if msg.strip() in http_errs:
        print('Tried getting data from KEGG, got HTTP err ' + msg)
        return True
    else: return False

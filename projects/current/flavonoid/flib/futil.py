import json
import os
import re

SEP = os.sep  # get the right slash. / for linux, \ for windows
JKEY = 'obj'

def get_json_data(file_name, key=None):
    """
    This function uses the python JSON library in order to parse JSON files into usable python objects. Can return
    lists or dictionaries, depending on the JSON file's structure.
    """
    if key is None: key = JKEY
    data = ''
    try:
        with open(file_name) as jsonFile:
            data = json.load(jsonFile)
    except FileNotFoundError:
        file_name = '..' + SEP + file_name
        with open(file_name) as jsonFile:
            data = json.load(jsonFile)
    return data[key]

def remove_dupes(dupe_list):
    """
    removes duplicate elements
    """
    unique_list = []  # creates an empty list
    for item in dupe_list:
        # adds item to empty list if it's not already in the list
        if item not in unique_list: unique_list.append(item)
    return unique_list

def unique_element_list(list_name, index):
    """
    find unique EC numbers so have a generic function and run it
    be careful as there are cases of one less item - use "last" to fix that problem here
    """
    original_index = index
    element_list = []
    for i in list_name:  # assigns the string "last" to the very last list in the list of lists
        if original_index == 'last': index = len(i) - 1  # finds unique EC num not in the list & adds it to the list
        if i[int(index)] not in element_list: element_list.append(i[int(index)])
    return element_list

def list_partition(seq, num):
    """
    This function takes in a list and then splits it into as many parts as specified by parameter num
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def write_readme(main_dir, readme, init_time, fasta_path, gene_path):
    """
    Creates ReadMe file
    """
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
    """
    This function takes a list or list of lists and then writes its contents to a file.
    """
    os.chdir(current)

    writedoc = open(output_dir, 'w')  # open file to be written
    for line in lists_to_write:
        for item in line:
            item = str(item).replace('\n', '')  # removes the new lines in each list of list
            if item == '':
                writedoc.write('-')  # if the list in the list of list is empty writes a dash
            else:
                writedoc.write(item)  # write the entry in the list of lists to the file
            writedoc.write(', ')  # tab delineated; use "," for csv files
    writedoc.write('\n')
    writedoc.close()

def write_append(path, content, write_over=None, skip=None):
    """
    This function takes a file name and the contents to be written to a file. If the file doesn't exist, it is
    created then written to. If it does exist, then it is appended to. The optional arg write_over is used for when
    the file exists and the content needs to be written over.
    """
    try:
        file = open(path, 'x')
        file.close()
        file = open(path, 'w')
        file.write(content)
        file.close()
    except FileExistsError:
        if write_over:
            file = open(path, 'w')
            file.write(content)
            file.close()
        elif skip:
            pass
        else:
            file = open(path, 'a')
            file.write(content)
            file.close()


def is_http_error(msg):
    """
    Checks if a string is an HTTP error.
    """
    if str(msg).strip() in HTTP_ERRS:
        return False
    else:
        return True

def init_dirs(main_dir, gene, fasta, chem):
    """
    Initializes directories for keggv2.py
    """
    # replaced WindowsError with OSError for more general usage. try to make data directories and handle any errors
    try:
        os.mkdir(main_dir)
    except OSError:
        pass
    try:
        os.mkdir(gene)
    except OSError or FileExistsError:
        pass
    try:
        os.mkdir(fasta)
    except OSError or FileExistsError:
        pass
    try:
        os.mkdir(chem)
    except OSError or FileExistsError:
        pass

def quick_fetch(pattern, line):
    """
    Fetch one item from re.findall and return as a string.
    """
    out = ''
    try:
        out = re.findall(pattern, line)[0]
    except IndexError:
        out = ''
    return out

def mult_replace(line, pairs):
    """
    This function makes multiple string replacements.
    """
    for pair in pairs: line = line.replace(pair[0], pair[1])
    return line

def skin(line):
    """
    Simple function that removes ALL whitespace from a string
    """
    return ''.join(line.split()).strip()

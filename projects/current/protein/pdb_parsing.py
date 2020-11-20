from Bio.PDB import PDBIO, MMCIFParser
from pathlib import Path
from plib.types import *
from threading import Lock, Thread

import datetime
import multiprocessing
import os
import sys
import urllib.error
import urllib.request
import math

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))  # allows for imports from directories at the same level
from lib.jsondata import *
from lib.util import *

pdb_objects_list = []
pdb_objects = {}
pdb_entries = []
pdb_basic_info = ''
total_pdb_output = ''

init_time = datetime.datetime.now()  # record start time of program

# make sure that the output directories exist before continuing
Path(out_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)

thread_lim = multiprocessing.cpu_count() - 1  # determine number of usable threads for the program

def main():
    ans = True
    print('Enter your choice:')
    while ans:
        print(
            '\n1. Download & parse all files' +
            '\n2. Run calculations on single file' +
            '\n0. Exit'
        )

        ans = int(input('\nSelection: '))
        if ans == 1:
            run_parse()
        elif ans == 2:
            to_use = input('\nEnter PDB ID: ').strip().upper()
            this_url = PART_URL + to_use + '.pdb'
            calc(this_url, to_use)
            pass

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken: ' + str(total_time))

def run_parse():
    split_pdb_list = list_partition(pdb_id_list, thread_lim)
    threads = []

    for chunks in split_pdb_list:
        tmp_thread = Thread(target=run_pdb_chunks, args=(chunks,))
        tmp_thread.start()
        threads.append(tmp_thread)

    for thread in threads: thread.join()  # wait for all threads to finish execution

    with open(formatted_out, 'w+') as out_file:
        out_file.write(total_pdb_output)
        out_file.close()

    with open(formatted_basic, 'w+')as out_basic:
        out_basic.write(pdb_basic_info)
        out_basic.close()

def run_pdb_chunks(chunks):
    for pdb_id in chunks:
        tst_url = PART_URL + pdb_id + '.pdb'
        pdb_stuff(tst_url, pdb_dir + pdb_id + '.pdb', pdb_id)

def calc(pdb_url, pdb_id):
    global total_pdb_output, pdb_entries, pdb_objects, pdb_basic_info
    pdb_stuff(pdb_url, pdb_dir + pdb_id + '.pdb', pdb_id, skip_download=True)
    ans = True
    print('Enter your choice: ')
    while ans:
        print(
            '\n1. Find distance between two atoms' +
            '\n2. Find name of any H within 1.2 angstroms of any O'
            '\n0. Back to main menu'
        )
        ans = int(input('\nSelection: '))
        entry = get_entry(pdb_id)
        if ans == 1:
            name1 = input('Name for first atom: ').strip().upper()
            name2 = input('Name for second atom: ').strip().upper()
            rec1 = Record()
            rec2 = Record()
            found = 0

            for rec in entry.records:
                if rec.a_name == name1:
                    rec1 = rec
                    found += 1
                elif rec.a_name == name2:
                    rec2 = rec
                    found += 1

            if found != 2: print('err finding atoms. had ' + str(found) + ' matches')
            else: print_atom_distance(rec1, rec2)
        elif ans == 2:
            find_hydrogen(entry)
            pass

def get_entry(pdb_id):
    global pdb_entries
    for entry in pdb_entries:
        if pdb_id == entry.pdb_id: return entry
    print('!!!WARNING!!!Could not find entry for PDB ID {}, returning empty entry'.format(pdb_id))

def print_atom_distance(re1: Record, re2: Record):
    dist = get_atom_distance(re1, re2)
    print('distance between {} ({}, {}, {}) & {} ({}, {}, {}) is: {}'.format(re1.a_name, re1.x, re1.y, re1.z,
                                                                             re2.a_name, re2.x, re2.y, re2.z, dist))

def get_atom_distance(re1: Record, re2: Record):
    return math.sqrt(((re2.x - re1.x) ** 2) + ((re2.y - re1.y) ** 2) + ((re2.z - re1.z) ** 2))

def find_hydrogen(entry: Entry):
    o_list = []
    h_list = []
    close_h = []
    for rec in entry.records:
        if rec.elem == 'O': o_list.append(rec)
        elif rec.elem == 'H': h_list.append(rec)

    for o_item in o_list:
        for h_item in h_list:
            if get_atom_distance(o_item, h_item) <= 1.20 and h_item not in close_h:
                print_atom_distance(o_item, h_item)
                close_h.append(h_item)

    close_str = ''
    o_str = ''
    h_str = ''
    for item in o_list: o_str += item.a_name + ' '
    for item in h_list: h_str += item.a_name + ' '
    for item in close_h: close_str += item.a_name + ' '
    print('O atoms: {}, H atoms: {}, H atoms within 1.2 angstroms of an O: {}'.format(o_str, h_str, close_str))


def pdb_stuff(url, path, pdb_id, skip_download=False):
    global total_pdb_output, pdb_entries, pdb_objects, pdb_basic_info
    lock_total = Lock()
    lock_entry = Lock()
    lock_obj = Lock()
    lock_basic = Lock()
    if not skip_download:
        if not os.path.exists(path):
            print('***PDB file not found for ' + pdb_id + ', starting download')

            try: urllib.request.urlretrieve(url, path)
            except urllib.error.HTTPError or urllib.error.URLError as e:
                print("!!!HTTP/URL error, couldn't get pdb file " + pdb_id + '. Error: ' + e.reason +
                      '. Will look for .cif file.')
                cif = path.replace('.pdb', '.cif')
                url = url.replace('.pdb', '.cif')
                cif_stuff(url, cif, path)

    try: file = open(path, 'r')
    except FileNotFoundError:
        print('!!!PDB file not found')
        return
    f_lines = file.readlines()
    tmp_entry = new_entry(lines=f_lines)
    tmp_entry.pdb_id = pdb_id

    for line in f_lines:
        if K_ATM in line or K_HAT in line:
            tmp_record = new_record(line=line, name=pdb_id)
            if tmp_record.lig_code.strip() in ligand_codes:
                with lock_total: total_pdb_output += tmp_record.important_str()
                with lock_obj:
                    if pdb_id in pdb_objects.keys(): pdb_objects[pdb_id].append(tmp_record)
                    else: pdb_objects[pdb_id] = [tmp_record]
                tmp_entry.records.append(tmp_record)
                # print(tmp_record.important_str())

    with lock_entry: pdb_entries.append(tmp_entry)
    simple_entry_print(pdb_id, tmp_entry.group, tmp_entry.ec_str, len(tmp_entry.records), tmp_entry.org_sci,
                       tmp_entry.org_taxid, tmp_entry.ex_sys)
    with lock_basic:
        pdb_basic_info += simple_entry_str(pdb_id, tmp_entry.group, tmp_entry.ec_str, len(tmp_entry.records),
                                           tmp_entry.org_sci, tmp_entry.org_taxid, tmp_entry.ex_sys) + '\n'

def cif_stuff(url, cif_path, pdb_path):
    try: urllib.request.urlretrieve(url, cif_path)
    except urllib.error.HTTPError or urllib.error.URLError as e:
        print("!!!HTTP or URL error, couldn't get " + url + '. Got error: ' + e.reason)
        return
    p = MMCIFParser()
    struc = p.get_structure('', cif_path)
    io = PDBIO()
    io.set_structure(struc)
    io.save(pdb_path)
    print('^^^SUCCESSFULLY CONVERTED CIF TO PDB')

def print_4v4d_pyg_chain_a(pdb_id, ligand, chain, line):
    if pdb_id == '4V4D' and ligand == 'PYG' and chain == 'A':
        print('[4V4D-PYG-Chain A] ' + line)

def simple_entry_print(pdb_id, group, ec_nums, num_records, org_name, tax_id, ex_sys):
    print(simple_entry_str(pdb_id, group, ec_nums, num_records, org_name, tax_id, ex_sys))

def simple_entry_str(pdb_id, group, ec_nums, num_records, org_name, tax_id, ex_sys):
    return pdb_id + ' || ' + group + ' || ' + org_name + ' || ' + tax_id + ' || ' + ex_sys + ' || ' + ec_nums + \
           ' || ' + str(num_records) + ' records'

if __name__ == '__main__':
    main()

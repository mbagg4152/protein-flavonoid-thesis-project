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

# make sure that the output directories exist before continuing
Path(out_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)

thread_lim = multiprocessing.cpu_count() - 1  # determine number of usable threads for the program

def main():
    # print(str(struct_planes))
    ans = True
    print('Enter your choice:')
    while ans:
        print('\n1. Download & parse all files'
              '\n2. Run calculations on single file'
              '\n0. Exit'
              )

        ans = int(input('\nSelection: '))
        if ans == 1:
            run_parse()
        elif ans == 2:
            to_use = input('\nEnter PDB ID: ').strip().upper()
            calc(to_use)
            pass


def run_parse():
    init_time = datetime.datetime.now()  # record start time of downloading and parsing
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
    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken to download & parse files: {}'.format(total_time))

def run_pdb_chunks(chunks):
    for pdb_id in chunks:
        tst_url = PART_URL + pdb_id + '.pdb'
        get_parse_pdbs(tst_url, pdb_dir + pdb_id + '.pdb', pdb_id)

def calc(pdb_id):
    pdb_url = PART_URL + pdb_id + '.pdb'
    global total_pdb_output, pdb_entries, pdb_objects, pdb_basic_info
    get_parse_pdbs(pdb_url, pdb_dir + pdb_id + '.pdb', pdb_id, skip_download=True)
    ans = True
    while ans:
        print('Enter your choice: '
              '\n1. Find distance between two atoms'
              '\n2. Find name of any H within 1.2 angstroms of any O'
              '\n3. Add new plane'
              '\n9. Change structure ID'
              '\n0. Back to main menu'
              )
        ans = int(input('\nSelection: '))
        struct = get_struct(pdb_id)
        if ans == 1:
            name1 = input('Name for first atom: ').strip().upper()
            name2 = input('Name for second atom: ').strip().upper()
            rec1 = Record()
            rec2 = Record()
            found = 0

            for rec in struct.records:
                if rec.atom == name1:
                    rec1 = rec
                    found += 1
                elif rec.atom == name2:
                    rec2 = rec
                    found += 1

            if found != 2: print('err finding atoms. had ' + str(found) + ' matches')
            else: print_atom_distance(rec1, rec2)
        elif ans == 2:
            find_hydrogen(struct)
            pass
        elif ans == 3:
            if pdb_id not in struct_planes: print('No plane definitions found for structure!')
            else: make_plane(struct)
        elif ans == 9:
            pdb_id = input('\nEnter PDB ID: ').strip().upper()
            pdb_url = PART_URL + pdb_id + '.pdb'
            get_parse_pdbs(pdb_url, pdb_dir + pdb_id + '.pdb', pdb_id, skip_download=True)
            struct = get_struct(pdb_id)

def get_struct(pdb_id):
    global pdb_entries
    for entry in pdb_entries:
        if pdb_id == entry.pdb_id: return entry
    print('!!!WARNING!!!Could not find entry for PDB ID {}, returning empty entry'.format(pdb_id))
    return Struct()

def print_atom_distance(re1: Record, re2: Record):
    dist = get_atom_distance(re1, re2)
    print('distance between {} ({}, {}, {}) & {} ({}, {}, {}) is: {}'.format(re1.atom, re1.x, re1.y, re1.z,
                                                                             re2.atom, re2.x, re2.y, re2.z, dist))

def get_atom_distance(re1: Record, re2: Record):
    try: return math.sqrt(((re2.x - re1.x) ** 2) + ((re2.y - re1.y) ** 2) + ((re2.z - re1.z) ** 2))
    except TypeError: return 'BAD VALUE!'

def find_hydrogen(entry: Struct):
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
    for item in o_list: o_str += item.atom + ' '
    for item in h_list: h_str += item.atom + ' '
    for item in close_h: close_str += item.atom + ' '
    print('O atoms: {}, H atoms: {}, H atoms within 1.2 angstroms of an O: {}'.format(o_str, h_str, close_str))

def make_plane(struct):
    for plane in struct.planes: print(plane.string())


def get_parse_pdbs(url, path, pdb_id, skip_download=False):
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
                get_convert_cifs(url, cif, path)

    try: file = open(path, 'r')
    except FileNotFoundError:
        print('!!!PDB file not found')
        return
    f_lines = file.readlines()
    tmp_entry = new_struct(lines=f_lines, pdb_id=pdb_id)

    for line in f_lines:
        if K_ATM in line or K_HAT in line:
            tmp_record = new_record(line=line, name=pdb_id)
            if tmp_record.lig_code.strip() in ligand_codes:
                with lock_total: total_pdb_output += tmp_record.info_str()
                with lock_obj:
                    if pdb_id in pdb_objects.keys(): pdb_objects[pdb_id].append(tmp_record)
                    else: pdb_objects[pdb_id] = [tmp_record]
                tmp_entry.records.append(tmp_record)

    for key in struct_planes:
        if key == pdb_id:
            for atoms in struct_planes[key]: tmp_entry.new_plane(atoms)
    with lock_entry: pdb_entries.append(tmp_entry)
    with lock_basic:
        pdb_basic_info += entry_str(pdb_id, tmp_entry.group, tmp_entry.ec_str, len(tmp_entry.records),
                                    tmp_entry.org_sci, tmp_entry.org_taxid, tmp_entry.ex_sys) + '\n'

def get_convert_cifs(url, cif_path, pdb_path):
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

def show_entry(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys):
    print(entry_str(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys))

def entry_str(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys):
    return '{} || {} || {} || {} || {} || {} || {} records'.format(pdb_id, grp, org, tax, ex_sys, ec_nums, num_rec)

if __name__ == '__main__':
    main()

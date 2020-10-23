from Bio.PDB import PDBIO, MMCIFParser
from plib.Types import *
from pathlib import Path
from threading import Lock, Thread
import datetime
import multiprocessing
import numpy
import sys
import urllib.request, urllib.error
import os
import time

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))  # allows for imports from directories at the same level

from lib.jsondata import *


pdb_objects_list = []
pdb_objects = {}
pdb_entries = []
total_pdb_output = ''

init_time = datetime.datetime.now()

Path(out_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(sasa_dir).mkdir(parents=True, exist_ok=True)

thread_lim = multiprocessing.cpu_count()


def main():
    chunked_pdb = numpy.array_split(numpy.array(pdb_id_list), thread_lim)
    threads = []
    print(str(len(chunked_pdb)))
    for chunks in chunked_pdb:
        tmp_thread = Thread(target=run_pdb_chunks, args=(chunks.tolist(),))
        tmp_thread.start()
        threads.append(tmp_thread)

    for thread in threads: thread.join()

    with open(formatted_out, 'w+') as out_file:
        out_file.write(total_pdb_output)
        out_file.close()

    for pdb_id in pdb_id_list:
        run_sasa(pdb_id, pdb_dir + pdb_id + '.pdb')

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken: ' + str(total_time))

def run_pdb_chunks(chunks):
    for pdb_id in chunks:
        tst_url = PART_URL + pdb_id + '.pdb'
        pdb_stuff(tst_url, pdb_dir + pdb_id + '.pdb', pdb_id)

def pdb_stuff(url, path, pdb_id):
    global total_pdb_output, pdb_entries, pdb_objects
    lock_total = Lock()
    lock_entry = Lock()
    lock_obj = Lock()
    if not os.path.exists(path):
        print('>>>PDB file not found for ' + pdb_id + ', starting download')

        try: urllib.request.urlretrieve(url, path)
        except urllib.error.HTTPError or urllib.error.URLError as e:
            print("<<<HTTP/URL error, couldn't pdb file " + pdb_id + '. Error: ' + e.reason +
                  '. Will look for .cif file.')
            cif = path.replace('.pdb', '.cif')
            url = url.replace('.pdb', '.cif')
            cif_stuff(url, cif, path)

    try: file = open(path, 'r')
    except FileNotFoundError:
        print('PDB file not found')
        return
    f_lines = file.readlines()
    tmp_entry = new_entry(lines=f_lines)

    for line in f_lines:
        if K_ATM in line or K_HAT in line:
            tmp_record = new_record(line, pdb_id)
            if tmp_record.ligand_code.strip() in ligand_codes:
                with lock_total: total_pdb_output += tmp_record.important_str()
                with lock_obj:
                    if pdb_id in pdb_objects.keys(): pdb_objects[pdb_id].append(tmp_record)
                    else: pdb_objects[pdb_id] = [tmp_record]

                tmp_entry.records.append(tmp_record)
                # print_4v4d_pyg_chain_a(pdb_id, tmp_record.ligand_code, tmp_record.chain_id, line)
    simple_entry_print(pdb_id, tmp_entry.group, tmp_entry.ec_nums, len(tmp_entry.records))
    with lock_entry: pdb_entries.append(tmp_entry)
    # run_sasa(pdb_id, path)

def cif_stuff(url, cif_path, pdb_path):
    try: urllib.request.urlretrieve(url, cif_path)
    except urllib.error.HTTPError or urllib.error.URLError as e:
        print("<<<HTTP or URL error, couldn't get " + url + '. Got error: ' + e.reason)
        return
    p = MMCIFParser()
    struc = p.get_structure('', cif_path)
    io = PDBIO()
    io.set_structure(struc)
    io.save(pdb_path)
    print('@@@SUCCESSFULLY CONVERTED CIF TO PDB')

def print_4v4d_pyg_chain_a(pdb_id, ligand, chain, line):
    if pdb_id == '4V4D' and ligand == 'PYG' and chain == 'A':
        print('[4V4D-PYG-Chain A] ' + line)

def simple_entry_print(pdb_id, group, ec_nums, num_records):
    print(pdb_id + ' | Class: ' + group + ' | EC: ' + str(ec_nums) + ' | No. Records: ' + str(num_records))

def run_sasa(pdb_id, path):
    try:
        os.mkdir(sasa_dir + pdb_id)
        time.sleep(1)
        cmd = 'cd; ' + \
              'cd ' + sasa_dir + pdb_id + ';' + \
              sasa + ' -m 4 -i ' + path + '> out.txt'
        os.system(cmd)
    except FileExistsError: pass
if __name__ == '__main__':
    main()

import os
import sys
import urllib.request, urllib.error
from pathlib import Path
from StringsAndConsts import *
from threading import Lock, Thread
import numpy
import datetime
import multiprocessing

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))
from Types import *
from lib.jsondata import *


pdb_objects_list = []
pdb_objects = {}
pdb_entries = []
total_pdb_output = ''

init_time = datetime.datetime.now()
print(str(init_time))

Path(out_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)

thread_lim = multiprocessing.cpu_count()
print(str(thread_lim))


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

    end_time = datetime.datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken: ' + str(total_time))


# def run_threads():
#     chunked_pdb = numpy.array_split(numpy.array(pdb_id_list), thread_lim)
#     threads = []
#
#     for chunks in chunked_pdb:
#         tmp_thread = Thread(target=run_pdb_chunks, args=(chunks.tolist(),))
#         tmp_thread.start()
#         threads.append(tmp_thread)
#
#     for thread in threads: thread.join()
#     return


# def update_thread_list(threads):
#     new_threads = [t for t in threads if t.is_alive()]
#     return new_threads
#
# def run_pdb(pdb_id):
#     tst_url = PART_URL + pdb_id + '.pdb'
#     pdb_stuff(tst_url, pdb_dir + pdb_id + '.pdb', pdb_id)


def run_pdb_chunks(chunks):
    for pdb_id in chunks:
        tst_url = PART_URL + pdb_id + '.pdb'
        pdb_stuff(tst_url, pdb_dir + pdb_id + '.pdb', pdb_id)

def pdb_stuff(url, path, pdb_id):
    global total_pdb_output, pdb_entries, pdb_objects
    lock_total = Lock()
    lock_entry = Lock()
    lock_obj = Lock()
    lock_print = Lock()
    if not os.path.exists(path):
        print('>>>PDB file not found for ' + pdb_id + ', starting download')

        try: urllib.request.urlretrieve(url, path)
        except urllib.error.HTTPError or urllib.error.URLError as e:
            print("<<<HTTP or URL error, couldn't get file for " + pdb_id + '. Got error: ' + e.reason +
                  '. Will look for .cif file.')
            # cif_stuff(url, path, pdb_id)
            return

    file = open(path, 'r')
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

    print(tmp_entry.pdb_id + ' | Class: ' + tmp_entry.group + ' | EC: ' + str(tmp_entry.ec_nums) +
          ' | No. Records: ' + str(len(tmp_entry.records)))

    with lock_entry: pdb_entries.append(tmp_entry)

def cif_stuff(url, path, pdb_id):
    print('')

if __name__ == '__main__':
    main()

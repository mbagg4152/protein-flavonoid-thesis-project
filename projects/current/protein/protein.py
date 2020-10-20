import os
import sys
import urllib.request
from pathlib import Path

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))
from Types import *
from lib.jsondata import *

partial_url = "https://files.rcsb.org/view/"
pdb_objects_list = []
pdb_objects = {}
pdb_entries = []
total_pdb_output = ''
cwd = os.getcwd()
out_dir = cwd + SEP + 'output' + SEP
formatted_out = out_dir + 'format_out_all.txt'
formatted_out_small = out_dir + 'format_out_partial.txt'
pdb_dir = out_dir + 'pdb_files' + SEP

Path(out_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)


def main():
    for pdb_id in pdb_id_list_short:
        print('looking at: ' + str(pdb_id))
        tst_url = partial_url + pdb_id + '.pdb'
        pdb_stuff(tst_url, pdb_dir + pdb_id + '.pdb', pdb_id)
    with open(formatted_out_small, 'w+') as out_file:
        out_file.write(total_pdb_output)
        out_file.close()


def pdb_stuff(url, path, pdb_id):
    global total_pdb_output
    urllib.request.urlretrieve(url, path)
    file = open(path, 'r')
    f_lines = file.readlines()
    tmp_entry = new_entry(lines=f_lines)

    for line in f_lines:
        if 'ATOM' in line or 'HETATM' in line:
            tmp_record = new_record(line, pdb_id)
            if tmp_record.ligand_code.strip() in ligand_codes:
                total_pdb_output += tmp_record.important_str()
                if pdb_id in pdb_objects.keys(): pdb_objects[pdb_id].append(tmp_record)
                else: pdb_objects[pdb_id] = [tmp_record]
            tmp_entry.records.append(tmp_record)
    print('tmp entry | PDB id: ' + tmp_entry.pdb_id + ' | header: ' + tmp_entry.header + ' | title: ' +
          tmp_entry.title + ' | ec nums: ' + str(tmp_entry.ec_nums) + ' | no. of (relevant) atom records: ' +
          str(len(tmp_entry.records)))
    pdb_entries.append(tmp_entry)


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

if __name__ == '__main__':
    main()

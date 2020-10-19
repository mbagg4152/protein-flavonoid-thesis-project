import os
import sys
import urllib.request
from Record import *

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))
from lib.jsondata import *


partial_url = "https://files.rcsb.org/view/"


def main():
    tst_url = partial_url + '6DHL.pdb'
    pdb_stuff(tst_url, '6DHL.pdb')


def pdb_stuff(url, name):
    urllib.request.urlretrieve(url, name)
    file = open(name, 'r')
    f_lines = file.readlines()
    out_str = ''
    for line in f_lines:
        if 'ATOM' in line or 'HETATM' in line:
            tmp_record = new_record(line)
            if tmp_record.lig_code.strip() in ligand_codes:
                out_str += tmp_record.relevant_str()
    print(out_str)

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

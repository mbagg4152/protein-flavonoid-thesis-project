import sys
import os
from pathlib import Path
import urllib.request, urllib.error

sep = os.sep
DS_PATH = os.getcwd() + sep + 'dr_sasa' + sep + 'dr_sasa_n' + sep + 'build' + sep + 'dr_sasa'
DS_OUT = '..' + sep + '..' + sep + '..' + sep + 'sasa_out'
DSO_REL = os.getcwd() + sep + 'sasa_out'
sys.path.append(os.getcwd().replace(sep + 'protein', ''))  # allows for imports from directories at the same level
from plib.strings_consts import *
from lib.jsondata import *

Path(DSO_REL).mkdir(parents=True, exist_ok=True)


def main():
    for pid in pdb_id_single:
        tst_url = PART_URL + pid + '.pdb'
        tmp_path = DSO_REL + sep + pid
        tmp_file = tmp_path + sep + pid + '.pdb'
        Path(tmp_path).mkdir(parents=True, exist_ok=True)
        fetch_run(tst_url, tmp_path, tmp_file)


def fetch_run(url, dir_path, file_path):
    cmd = 'cd ' + dir_path + ';' + \
          'wget -O ' + file_path + ' ' + url + ';' + \
          DS_PATH + ' -m 4 -i ' + file_path
    os.system(cmd)


if __name__ == '__main__':
    main()

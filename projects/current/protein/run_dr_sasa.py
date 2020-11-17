import sys
import os
from pathlib import Path
import urllib.request, urllib.error

sep = os.sep
DS_PATH = os.getcwd() + sep + 'dr_sasa' + sep + 'dr_sasa_n' + sep + 'build' + sep + 'dr_sasa'
DS_OUT_M4 = '..' + sep + '..' + sep + '..' + sep + 'sasa_out_mode_4'
DS_OUT_M0 = '..' + sep + '..' + sep + '..' + sep + 'sasa_out_mode_0'
DSO_REL_M4 = os.getcwd() + sep + 'sasa_out_mode_4'
DSO_REL_M0 = os.getcwd() + sep + 'sasa_out_mode_0'
sys.path.append(os.getcwd().replace(sep + 'protein', ''))  # allows for imports from directories at the same level
from plib.strings_consts import *
from lib.jsondata import *

Path(DSO_REL_M4).mkdir(parents=True, exist_ok=True)


def main():
    for pid in pdb_id_single:
        tst_url = PART_URL + pid + '.pdb'
        tmp_path = DSO_REL_M4 + sep + pid
        tmp_path2 = DSO_REL_M0 + sep + pid
        tmp_file = tmp_path + sep + pid + '.pdb'
        tmp_file2 = tmp_path2 + sep + pid + '.pdb'
        Path(tmp_path).mkdir(parents=True, exist_ok=True)
        fetch_run(tst_url, tmp_path, tmp_path2, tmp_file, tmp_file2)


def fetch_run(url, dir_path1, dir_path2, file_path1, file_path2):
    cmd1 = 'cd ' + dir_path1 + ';' + \
           'wget -O ' + file_path1 + ' ' + url + ';' + \
           DS_PATH + ' -m 4 -i ' + file_path1
    cmd2 = 'cd ' + dir_path2 + ';' + \
           'wget -O ' + file_path2 + ' ' + url + ';' + \
           DS_PATH + ' -m 4 -i ' + file_path2
    os.system(cmd1)


if __name__ == '__main__':
    main()

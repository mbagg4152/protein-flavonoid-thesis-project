import sys
import os
from os.path import exists
from pathlib import Path
import platform
import urllib.request as url_req
import urllib.error as url_err

sep = os.sep
DS_DIR = os.getcwd() + sep + 'dr_sasa'
DS_DIR_N = DS_DIR + sep + 'dr_sasa_n'
BUILD_DIR = DS_DIR_N + sep + 'build'
DS_PATH = BUILD_DIR + sep + 'dr_sasa'
DS_OUT_M4 = '..' + sep + '..' + sep + '..' + sep + 'm4'
DS_OUT_M0 = '..' + sep + '..' + sep + '..' + sep + 'm0'
DSO_REL_M4 = os.getcwd() + sep + 'm4'
DSO_REL_M0 = os.getcwd() + sep + 'm0'
sys.path.append(os.getcwd().replace(sep + 'protein', ''))  # allows for imports from directories at the same level
from plib.strings_consts import *

try: from json_objects import *
except ImportError:
    # allows for imports from directories at the same level
    sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))
    from sharedlib.json_objects import *

Path(DSO_REL_M4).mkdir(parents=True, exist_ok=True)
Path(DSO_REL_M0).mkdir(parents=True, exist_ok=True)

def main():
    check_install()
    if not exists(DS_DIR) or not exists(DS_DIR_N) or not exists(BUILD_DIR) or not exists(DS_PATH):
        print('Something went wrong in setup')
        return

    for pid in pdb_id_list_short:
        tst_url = PART_URL + pid + '.pdb'
        tmp_path = DSO_REL_M4 + sep + pid
        tmp_path2 = DSO_REL_M0 + sep + pid
        tmp_file = tmp_path + sep + pid + '.pdb'
        tmp_file2 = tmp_path2 + sep + pid + '.pdb'
        Path(tmp_path).mkdir(parents=True, exist_ok=True)
        Path(tmp_path2).mkdir(parents=True, exist_ok=True)
        fetch_run(tst_url, tmp_path, tmp_path2, tmp_file, tmp_file2)

def check_install():
    if not exists(DS_DIR) or not exists(DS_DIR_N) or not exists(BUILD_DIR) or not exists(DS_PATH):
        if platform.system() != 'Linux':
            print('* This is a quick install only tested for Ubuntu (Linux). \n* Please make sure that for other '
                  'installations, that the dr_sasa executable/binary is in protein/dr_sasa/dr_sasa_n/build/\n'
                  '* Windows & Mac binaries can be found here: https://github.com/nioroso-x3/dr_sasa_n/releases\n'
                  '* Install instructions here: https://github.com/nioroso-x3/dr_sasa_n/blob/master/INSTALL\n'
                  '* dr_sasa home page: http://schuellerlab.org/dr_sasa/\n')
        else:
            cmd = 'mkdir dr_sasa; cd dr_sasa; git clone https://github.com/nioroso-x3/dr_sasa_n.git;' \
                  'cd dr_sasa_n/; mkdir build; cd build/; cmake ../; make; cd ../..;'
            try:
                import apt
            except ImportError:
                print('This machine does not have the apt package manager and cannot check for program installs\n'
                      'In case the install doesn\'t work, read the following:\n'
                      '* Make sure git, make and cmake are installed.\n'
                      '* Binaries found here: https://github.com/nioroso-x3/dr_sasa_n/releases\n'
                      '* Install instructions here: https://github.com/nioroso-x3/dr_sasa_n/blob/master/INSTALL\n'
                      '* Linux install instructions here: https://github.com/nioroso-x3/dr_sasa_n/blob/master/INSTALL')
                os.system(cmd)
                return
            cache = apt.Cache()
            missing = 0
            missing_str = ''
            if not cache['git'].is_installed:
                missing += 1
                missing_str += 'sudo apt-get -y install git\n'
            if not cache['cmake'].is_installed:
                missing += 1
                missing_str += 'sudo apt-get -y install cmake\n'
            if not cache['make'].is_installed:
                missing += 1
                missing_str += 'sudo apt-get -y install make\n'
            if missing > 0:
                print('Type the following command(s) into the terminal in order for the quick setup to work:\n' +
                      missing_str + 'Note: if you are using another linux distribution replace apt-get with the '
                                    'appropriate package manager')
                return
            print('Need to download and setup dr sasa')

            os.system(cmd)

def fetch_run(url, dir_path_m4, dir_path_m0, file_path_m4, file_path_m0):
    cmd_m0_1 = 'cd ' + dir_path_m0
    cmd_m0_2 = cmd_m0_1 + ';' + DS_PATH + ' -m 0 -i ' + file_path_m0
    cmd_m4_1 = 'cd ' + dir_path_m4
    cmd_m4_2 = cmd_m4_1 + ';' + DS_PATH + ' -m 4 -i ' + file_path_m4

    os.system(cmd_m0_1)
    try:
        url_req.urlretrieve(url, file_path_m0)
        print('running mode 0')

        os.system(cmd_m0_2)
    except (url_err.HTTPError, url_err.URLError): print('Could\'t get file from ' + url + ', will not run dr sasa in mode 0 for this pdb')
    os.system(cmd_m4_1)
    try:
        url_req.urlretrieve(url, file_path_m4)
        print('running mode 4')

        os.system(cmd_m4_2)
    except (url_err.HTTPError, url_err.URLError): print('Could\'t get file from ' + url + ', will not run dr sasa in mode 4 for this pdb')

if __name__ == '__main__':
    main()

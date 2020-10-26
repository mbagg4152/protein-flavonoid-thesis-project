import datetime
import os
import re
import sys
import threading
import time
import warnings
import logging.config

sys.path.append(os.getcwd().replace(os.sep + 'flavonoid', ''))

from bioservices.kegg import KEGG
from lib.jsondata import *
from lib.datatypes import *
from lib.pathstrings import *
from lib.compoundinfo import *

init_time = datetime.datetime.now()
kegg = KEGG()

path_fasta = ' '
path_gene = ' '
path_chem = ' '
main_dir = ' '
cwd = ' '


def main():
    pass


def make_dirs():
    global cwd, main_dir, path_chem, path_fasta, path_gene
    decision = ''
    if len(sys.argv) > 1: decision = sys.argv[1]
    else:
        print("No directory name supplied in args, defaulting to directory 'data'. "
              "Supply directory name in terminal using 'python3 kegg-prog.py dir_name'")
        decision = 'data'
    # create sub dirs for gene, FASTA & chemical data
    cwd = os.getcwd() + SEP
    main_dir = cwd + decision
    path_chem = main_dir + CHEM_DIR
    path_fasta = main_dir + FASTA_DIR
    path_gene = main_dir + GENE_DIR

    # replaced WindowsError with OSError for more general usage. try to make data directories and handle any errors
    try: os.mkdir(main_dir)
    except OSError: pass
    try: os.mkdir(path_gene)
    except OSError or FileExistsError: pass
    try: os.mkdir(path_fasta)
    except OSError or FileExistsError: pass
    try: os.mkdir(path_chem)
    except OSError or FileExistsError: pass


if __name__ == '__main__':
    main()

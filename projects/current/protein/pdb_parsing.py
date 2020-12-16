from datetime import datetime
from math import acos, degrees
from pathlib import Path
from plib.types import *
from threading import Lock, Thread
from urllib import request
import urllib.error
import multiprocessing
from plib.pdb_util import *

try:
    import matplotlib.pyplot as plt
    from Bio.PDB import PDBIO, MMCIFParser
except ImportError:
    print('Program needs matplotlib and biopython to work. One or both are missing. In the terminal, try:\n'
          '`pip3 install -U matplotlib` or `pip install -U matplotlib`\n'
          '`pip3 install biopython` or `pip install biopython`')
    exit(1)
try:
    from json_objects import *
    from util import *
except ImportError:
    # allows for imports from directories at the same level
    sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))
    from lib.json_objects import *
    from lib.util import *
plt.rc('xtick', labelsize=4)  # set font size for ticks on x axis for pyplot
plt.rc('ytick', labelsize=4)  # set font size for ticks on y axis for pyplot
pdb_objects_list = []
pdb_entries = []
pdb_objects = {}
pdb_basic_info = ''
total_pdb_output = ''

MAKE_GRAPHS = False

# make sure that the output directories exist before continuing
Path(out_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(image_dir).mkdir(parents=True, exist_ok=True)
thread_lim = multiprocessing.cpu_count() - 1  # determine number of usable threads for the program

def main():
    run_file_calc_thread()
    # run_parse()

def run_file_calc_thread():
    test_pdb_ids = get_json_data(FN_LIG_TESTS)
    for pdb in test_pdb_ids: file_calculations(pdb)

def run_parse():
    """For each of the PDB IDs, the PDB file will be downloaded if missing, then each PDB file is parsed accordingly."""
    init_time = datetime.now()  # record start time of downloading and parsing
    split_pdb_list = list_partition(pdb_id_list, thread_lim)
    threads = []

    for sub_list in split_pdb_list:
        tmp_thread = Thread(target=run_on_thread, args=(sub_list,))
        tmp_thread.start()
        threads.append(tmp_thread)

    for thread in threads: thread.join()  # wait for all threads to finish execution

    with open(formatted_out, 'w+') as out_file:
        out_file.write(total_pdb_output)
        out_file.close()

    with open(formatted_basic, 'w+')as out_basic:
        out_basic.write(pdb_basic_info)
        out_basic.close()
    end_time = datetime.now()
    total_time = end_time - init_time
    print('\ntotal time taken to download & parse files: {}'.format(total_time))

def run_on_thread(items):
    for pdb_id in items:
        tst_url = PART_URL + pdb_id + '.pdb'
        get_parse_pdbs(tst_url, pdb_dir + pdb_id + '.pdb', pdb_id)

def file_calculations(pdb_id):
    this_url = PART_URL + pdb_id + '.pdb'
    global total_pdb_output, pdb_entries, pdb_objects, pdb_basic_info
    get_parse_pdbs(this_url, pdb_test_dir + pdb_id + '.pdb', pdb_id, skip_download=True)

def find_hydrogen(entry: Struct):
    """Find every H atom that is within 1.2 angstroms of any O atom."""
    oxygens, hydrogens, selected_hydrogens = [], [], []
    for rec in entry.records:
        if rec.elem == 'O': oxygens.append(rec)
        elif rec.elem == 'H': hydrogens.append(rec)

    for o_item in oxygens:
        for h_item in hydrogens:
            if atom_distance(o_item, h_item) <= 1.20 and h_item not in selected_hydrogens:
                show_atom_distance(o_item, h_item)
                selected_hydrogens.append(h_item)

    close_str, o_str, h_str = '', '', ''
    for item in oxygens: o_str += item.atom + ' '
    for item in hydrogens: h_str += item.atom + ' '
    for item in selected_hydrogens: close_str += item.atom + ' '
    print('O atoms: {}, H atoms: {}, H atoms within 1.2 angstroms of an O: {}'.format(o_str, h_str, close_str))

def plane_operations(struct):
    if len(struct.rings) < 2:
        print('There are not enough planes to complete this action')
        return

    print('Struct {} has the following planes: '.format(struct.pdb_id))
    plane_count, p_index0, p_index1 = 0, 0, 1
    for plane in struct.rings:
        print('Plane {}. {}'.format(plane_count, plane.string()))
        plane_count += 1

    if len(struct.rings) == 2: print('Struct has only two planes, these will be used in the calculations')
    else: p_index0, p_index1 = int(input('\nFirst plane: ')), int(input('Second plane: '))
    try: plane0, plane1 = struct.rings[p_index0], struct.rings[p_index1]
    except IndexError:
        print('Incorrect index value given')
        return

    Path(image_dir + SEP + struct.pdb_id).mkdir(parents=True, exist_ok=True)
    if MAKE_GRAPHS:
        quick_plot(plane0, image_dir + SEP + struct.pdb_id + SEP + 'plane' + str(p_index0) + '_chain_' + plane0.chain)
        quick_plot(plane1, image_dir + SEP + struct.pdb_id + SEP + 'plane' + str(p_index1) + '_chain_' + plane1.chain)
    eqn0, eqn1 = plane0.eqn, plane1.eqn
    print('Plane {}. {} || {}'.format(p_index0, eqn0.string(), eqn0.func_form()))
    print('Plane {}. {} || {}'.format(p_index1, eqn1.string(), eqn1.func_form()))
    angle = plane_angles(eqn0, eqn1)
    print('the angle between the planes is {:6.3f} degrees'.format(angle))

def plane_angles(e0: Equation, e1: Equation):
    top = abs((e0.a * e1.a) + (e0.b * e1.b) + (e0.c * e1.c))
    bottom = (sqrt((e0.a ** 2) + (e0.b ** 2) + (e0.c ** 2))) * (sqrt((e1.a ** 2) + (e1.b ** 2) + (e1.c ** 2)))
    return degrees(acos(top / bottom))

def get_parse_pdbs(url, path, pdb_id, skip_download=False):
    global total_pdb_output, pdb_entries, pdb_objects, pdb_basic_info
    lock_total, lock_entry, lock_obj, lock_basic = Lock(), Lock(), Lock(), Lock()

    if not skip_download:
        if not os.path.exists(path):
            print('***PDB file not found for ' + pdb_id + ', starting download')

            try:
                urllib.request.urlretrieve(url, path)
                print('✓✓✓PDB ' + pdb_id + ' downloaded.')
            except urllib.error.URLError or urllib.error.HTTPError:
                print("!!!HTTP/URL error, couldn't get pdb file " + pdb_id + '. Will look for .cif file.')
                cif = path.replace('.pdb', '.cif')
                url = url.replace('.pdb', '.cif')
                get_convert_cifs(url, cif, path)
        else: print('$$$PDB file for ' + pdb_id + ' exists')

    try: file = open(path, 'r')
    except FileNotFoundError:
        print('!!!PDB file not found')
        return
    f_lines = file.readlines()
    tmp_entry = new_struct(lines=f_lines, pdb_id=pdb_id)
    tmp_entry.ligands = []
    for line in f_lines:
        if K_ATM in line or K_HAT in line:
            rec = new_record(line=line, name=pdb_id)

            if skin(rec.lig_code) in ligand_codes:
                with lock_total: total_pdb_output += rec.info_str()
                with lock_obj:
                    if pdb_id in pdb_objects.keys(): pdb_objects[pdb_id].append(rec)
                    else: pdb_objects[pdb_id] = [rec]
                tmp_entry.records.append(rec)
                if skin(rec.lig_code) not in tmp_entry.ligands: tmp_entry.ligands.append(skin(rec.lig_code))

    for lig in struct_rings:

        if lig in tmp_entry.ligands:
            rings = struct_rings.get(lig)
            for ring in rings:
                good_recs = []
                for atom in rings.get(ring):
                    good_recs.extend([o for o in tmp_entry.records if o.atom == atom and o.pdb_id == pdb_id])
                if len(good_recs):
                    print('!!!!matched records from {} for ligand {} ring {}'.format(pdb_id, lig, ring))
                    for gr in good_recs:
                        gr.show()

    with lock_entry: pdb_entries.append(tmp_entry)
    with lock_basic:
        pdb_basic_info += entry_str(pdb_id, tmp_entry.group, tmp_entry.ec_str, len(tmp_entry.records),
                                    tmp_entry.org_sci, tmp_entry.org_taxid, tmp_entry.ex_sys) + '\n'

if __name__ == '__main__':
    main()

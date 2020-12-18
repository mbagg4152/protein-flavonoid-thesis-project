from datetime import datetime
from math import acos, degrees
from pathlib import Path
from plib.types import *
from threading import Lock, Thread
from urllib import request
import urllib.error
import multiprocessing
from plib.pdb_util import *
from itertools import groupby

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
    with open(formatted_out, 'w+') as out_file:
        out_file.write(total_pdb_output)
        out_file.close()

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

def find_hydrogen(entry: ProtStruct):
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
    t_struct = new_struct(lines=f_lines, pdb_id=pdb_id)
    t_struct.ligands = []

    for line in f_lines:
        if K_ATM in line or K_HAT in line:
            rec = new_record(line=line, name=pdb_id)

            if skin(rec.lig_code) in ligand_codes:
                with lock_total: total_pdb_output += rec.info_str()
                with lock_obj:
                    if pdb_id in pdb_objects.keys(): pdb_objects[pdb_id].append(rec)
                    else: pdb_objects[pdb_id] = [rec]
                if isinstance(rec, AtomRec): t_struct.records.append(rec)
                if skin(rec.lig_code) not in t_struct.ligands: t_struct.ligands.append(skin(rec.lig_code))
    for lig in t_struct.ligands:
        if lig in list(struct_rings.keys()):
            rings = struct_rings.get(lig)
            for ring in rings:
                good_recs = []
                for atom_name in rings.get(ring):
                    for recs in t_struct.records:
                        # print('lig ' + lig)
                        if isinstance(recs, AtomRec):
                            try:
                                if recs.atom == atom_name and recs.pdb_id == pdb_id and recs.lig_code == lig: good_recs.append(recs)
                            except AttributeError: print('record was type {}'.format(str(type(recs))))
                if len(good_recs):
                    seq_list = list(set([o.seq for o in good_recs]))
                    chain_list = list(set([o.chain for o in good_recs]))
                    print(str(chain_list))
                    for seq in seq_list:
                        passed_recs = []
                        for gr in good_recs:
                            if gr.seq == seq: passed_recs.append(gr)

                        if len(passed_recs):
                            for ch in chain_list:
                                sect = []
                                names = []
                                for pr in passed_recs:
                                    if pr.chain.strip() == ch.strip():
                                        sect.append(pr)
                                        names.append(pr.atom)

                                if len(sect):
                                    try: t_struct.new_ring(names, sect[0].lig_code, sect[0].chain, sect[0].seq, ring)
                                    except AttributeError: pass
        for rng in t_struct.rings:
            if not isinstance(rng, Ring): t_struct.rings.remove(rng)
            else: print(pdb_id + ' ' + rng.string())
    with lock_entry: pdb_entries.append(t_struct)

def get_struct(pdb_id):
    """Find structure from global list of processed PDB files using PDB ID."""
    for entry in pdb_entries:
        if pdb_id == entry.pdb_id: return entry
    print('!!!WARNING!!!Could not find entry for PDB ID {}, returning empty entry'.format(pdb_id))
    return ProtStruct()

if __name__ == '__main__':
    main()

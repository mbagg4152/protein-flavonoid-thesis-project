from Bio.PDB import PDBIO, MMCIFParser
from datetime import datetime
from math import acos, degrees
from pathlib import Path
from plib.types import *
from threading import Lock, Thread
from urllib import request
from urllib.error import HTTPError, URLError
import matplotlib.pyplot as plt
import multiprocessing

plt.rc('xtick', labelsize=4)  # set font size for ticks on x axis for pyplot
plt.rc('ytick', labelsize=4)  # set font size for ticks on y axis for pyplot

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))  # allows for imports from directories at the same level
from lib.jsondata import *
from lib.util import *

pdb_objects_list = []
pdb_entries = []
pdb_objects = {}
pdb_basic_info = ''
total_pdb_output = ''

MAKE_GRAPHS = True

# make sure that the output directories exist before continuing
Path(out_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(image_dir).mkdir(parents=True, exist_ok=True)
thread_lim = multiprocessing.cpu_count() - 1  # determine number of usable threads for the program

def main():
    ans = True
    while ans:
        print('\n1. Download & parse all files'
              '\n2. Run calculations on single file'
              '\n0. Exit'
              '\n\n---------------------------------------------------------------------------------\n')
        ans = int(input('\nChoice: '))
        if ans == 1: run_parse()
        elif ans == 2:
            # to_use = input('\nEnter PDB ID: ').strip().upper()
            to_use = '2NNL'
            file_calculations(to_use)
            pass

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
    pdb_url = PART_URL + pdb_id + '.pdb'
    global total_pdb_output, pdb_entries, pdb_objects, pdb_basic_info
    get_parse_pdbs(pdb_url, pdb_dir + pdb_id + '.pdb', pdb_id, skip_download=True)
    ans = True
    while ans:
        print('1. Find distance between two atoms'
              '\n2. Find name of any H within 1.2 angstroms of any O'
              '\n3. Find angle between two planes'
              '\n9. Change structure ID'
              '\n0. Back to main menu'
              '\n\n---------------------------------------------------------------------------------\n')

        ans = int(input('Choice: '))
        struct = get_struct(pdb_id)

        if ans == 1:  # Find distance between two atoms
            name1, name2 = input('First atom: ').strip().upper(), input('Second atom: ').strip().upper()
            rec1, rec2 = Record(), Record()
            found = 0
            for rec in struct.records:
                if rec.atom == name1: rec1, found = rec, found + 1
                elif rec.atom == name2: rec2, found = rec, found + 1
            if found != 2: print('err finding atoms. had ' + str(found) + ' matches')
            else: show_atom_distance(rec1, rec2)

        elif ans == 2: find_hydrogen(struct)  # Find name of any H within 1.2 angstroms of any O

        elif ans == 3:  # Find angle between two planes
            if pdb_id not in struct_planes: print('No plane definitions found for structure!')
            else: plane_operations(struct)

        elif ans == 9:  # Change structure ID
            pdb_id = input('\nEnter PDB ID: ').strip().upper()
            pdb_url = PART_URL + pdb_id + '.pdb'
            get_parse_pdbs(pdb_url, pdb_dir + pdb_id + '.pdb', pdb_id, skip_download=True)
            struct = get_struct(pdb_id)

        print('\n---------------------------------------------------------------------------------\n')

def get_struct(pdb_id):
    """Find structure from global list of processed PDB files using PDB ID."""
    for entry in pdb_entries:
        if pdb_id == entry.pdb_id: return entry
    print('!!!WARNING!!!Could not find entry for PDB ID {}, returning empty entry'.format(pdb_id))
    return Struct()

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
    if len(struct.planes) < 2:
        print('There are not enough planes to complete this action')
        return

    print('Struct {} has the following planes: '.format(struct.pdb_id))
    plane_count, p_index0, p_index1 = 0, 0, 1
    for plane in struct.planes:
        print('Plane {}. {}'.format(plane_count, plane.string()))
        plane_count += 1

    if len(struct.planes) == 2: print('Struct has only two planes, these will be used in the calculations')
    else: p_index0, p_index1 = int(input('\nFirst plane: ')), int(input('Second plane: '))
    try: plane0, plane1 = struct.planes[p_index0], struct.planes[p_index1]
    except IndexError:
        print('Incorrect index value given')
        return

    Path(image_dir + SEP + struct.pdb_id).mkdir(parents=True, exist_ok=True)
    if MAKE_GRAPHS:
        quick_plot(plane0, image_dir + SEP + struct.pdb_id + SEP + 'plane' + str(p_index0) + '_' + plane0.chain)
        quick_plot(plane1, image_dir + SEP + struct.pdb_id + SEP + 'plane' + str(p_index1) + '_' + plane1.chain)
    eqn0, eqn1 = plane0.eqn, plane1.eqn
    print('Plane {}. {} || {}'.format(p_index0, eqn0.string(), eqn0.func_form()))
    print('Plane {}. {} || {}'.format(p_index1, eqn1.string(), eqn1.func_form()))
    angle = plane_angles(eqn0, eqn1)
    print('the angle between the planes is {:6.3f} degrees'.format(angle))

def plane_angles(e0: Eqn, e1: Eqn):
    top = abs((e0.a * e1.a) + (e0.b * e1.b) + (e0.c * e1.c))
    bottom = (sqrt((e0.a ** 2) + (e0.b ** 2) + (e0.c ** 2))) * (sqrt((e1.a ** 2) + (e1.b ** 2) + (e1.c ** 2)))
    return degrees(acos(top / bottom))

def same_atom(atom1, atom2):
    return (atom1.x == atom2.x) and (atom1.y == atom2.y) and (atom1.z == atom2.z)

def quick_plot(plane, path):
    """Output figures which show the planes and the atoms in a graph at different viewing angles."""
    atoms = plane.atoms
    atoms.sort(key=lambda a: a.from_center)

    # make list of x, y & z coordinates using each of the atoms
    x, y, z = list((atm.x for atm in atoms)), list((atm.y for atm in atoms)), list((atm.z for atm in atoms))
    fig = plt.figure(figsize=(10, 10))  # create base matplotlib figure

    # get evenly spaced numbers over a specified interval for x & y
    x_line, y_line = np.linspace(min(x), max(x), int(max(z)) + 1), np.linspace(min(y), max(y), int(max(z)) + 1)
    xs, ys = np.meshgrid(x_line, y_line)  # get coordinate matrices from coordinate vectors for x & y values
    vals = plane.eqn.func_form_tup()  # convert to function format
    zs = vals[0] * xs + vals[1] * ys + vals[2]  # use values from function format to get z

    # rotations & subplot grid placements
    angles, plots = [0, 45, 135, 0, 45, 135, 0, 45, 135], [331, 332, 333, 334, 335, 336, 337, 338, 339]
    views = [20, 20, 20, 60, 60, 60, 5, 5, 5]  # graph camera angle
    for i in range(0, len(plots)):
        ax = fig.add_subplot(plots[i], projection='3d')  # add subplot to figure
        ax.scatter(x, y, z, c='r')  # graph the atoms
        ax.plot_surface(xs, ys, zs, alpha=0.2)  # make it so that the points aren't covered by the plane
        ax.view_init(views[i], angles[i])  # make subplot at specific view & angle
    fig.savefig(path + '.png', dpi=200)

def get_parse_pdbs(url, path, pdb_id, skip_download=False):
    global total_pdb_output, pdb_entries, pdb_objects, pdb_basic_info
    lock_total, lock_entry, lock_obj, lock_basic = Lock(), Lock(), Lock(), Lock()

    if not skip_download:
        if not os.path.exists(path):
            print('***PDB file not found for ' + pdb_id + ', starting download')

            try: request.urlretrieve(url, path)
            except HTTPError or URLError as e:
                print("!!!HTTP/URL error, couldn't get pdb file " + pdb_id + '. Error: ' + e.reason +
                      '. Will look for .cif file.')
                cif = path.replace('.pdb', '.cif')
                url = url.replace('.pdb', '.cif')
                get_convert_cifs(url, cif, path)

    try: file = open(path, 'r')
    except FileNotFoundError:
        print('!!!PDB file not found')
        return
    f_lines = file.readlines()
    tmp_entry = new_struct(lines=f_lines, pdb_id=pdb_id)

    for line in f_lines:
        if K_ATM in line or K_HAT in line:
            tmp_record = new_record(line=line, name=pdb_id)
            if tmp_record.lig_code.strip() in ligand_codes:
                with lock_total: total_pdb_output += tmp_record.info_str()
                with lock_obj:
                    if pdb_id in pdb_objects.keys(): pdb_objects[pdb_id].append(tmp_record)
                    else: pdb_objects[pdb_id] = [tmp_record]
                tmp_entry.records.append(tmp_record)

    for key in struct_planes:
        if key == pdb_id:
            for lig in struct_planes[key]:
                for chain in struct_planes[key][lig]:
                    for atoms in struct_planes[key][lig][chain]:
                        print(str(atoms))
                        tmp_entry.new_plane(atoms, lig, chain)

    with lock_entry: pdb_entries.append(tmp_entry)
    with lock_basic:
        pdb_basic_info += entry_str(pdb_id, tmp_entry.group, tmp_entry.ec_str, len(tmp_entry.records),
                                    tmp_entry.org_sci, tmp_entry.org_taxid, tmp_entry.ex_sys) + '\n'

def get_convert_cifs(url, cif_path, pdb_path):
    try: request.urlretrieve(url, cif_path)
    except HTTPError or URLError as e:
        print("!!!HTTP or URL error, couldn't get " + url + '. Got error: ' + e.reason)
        return
    p = MMCIFParser()
    struc = p.get_structure('', cif_path)
    io = PDBIO()
    io.set_structure(struc)
    io.save(pdb_path)
    print('^^^SUCCESSFULLY CONVERTED CIF TO PDB')

def show_entry(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys):
    print(entry_str(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys))

def entry_str(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys):
    return '{} || {} || {} || {} || {} || {} || {} records'.format(pdb_id, grp, org, tax, ex_sys, ec_nums, num_rec)

if __name__ == '__main__':
    main()

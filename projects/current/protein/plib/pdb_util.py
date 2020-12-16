from types import *

def get_struct(pdb_id):
    """Find structure from global list of processed PDB files using PDB ID."""
    for entry in pdb_entries:
        if pdb_id == entry.pdb_id: return entry
    print('!!!WARNING!!!Could not find entry for PDB ID {}, returning empty entry'.format(pdb_id))
    return Struct()

def same_atom(atom1, atom2):
    return (atom1.x == atom2.x) and (atom1.y == atom2.y) and (atom1.z == atom2.z)


def quick_plot(plane, path):
    """Output figures which show the planes and the atoms in a graph at different viewing angles."""
    print('Making graphs....')
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

def get_convert_cifs(url, cif_path, pdb_path):
    try: request.urlretrieve(url, cif_path)
    except urllib.error.URLError or urllib.error.HTTPError as e:
        print("!!!HTTP or URL error, couldn't get " + url + '.')
        return
    try:
        p = MMCIFParser()
        struc = p.get_structure('', cif_path)
        io = PDBIO()
        io.set_structure(struc)
        io.save(pdb_path)
        print('^^^SUCCESSFULLY CONVERTED CIF TO PDB')
    except TypeError: print('Problem making pdb file')

def show_entry(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys):
    print(entry_str(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys))

def entry_str(pdb_id, grp, ec_nums, num_rec, org, tax, ex_sys):
    return '{} || {} || {} || {} || {} || {} || {} records'.format(pdb_id, grp, org, tax, ex_sys, ec_nums, num_rec)

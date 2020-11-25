from itertools import combinations
from math import sqrt
from plib.strings_consts import *
import sys

try:
    from scipy.optimize import leastsq
    import numpy as np
except ImportError:
    print('Program needs scipy and numpy to work. One or both are missing. In the terminal, try:\n'
          '`pip3 install scipy` or `pip install scipy`\n'
          '`pip3 install numpy` or `pip install numpy`')
    exit(1)

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))  # allows for imports from directories at the same level
from lib.util import *

class Atom:
    def __init__(self, name=DS, elem=DS, x=DF, y=DF, z=DF, from_center=DF, ligand=DS, chain=DS):
        self.name = name
        self.elem = elem
        self.x = x
        self.y = y
        self.z = z
        self.from_center = from_center
        self.ligand = ligand
        self.chain = chain
        self.from_origin()

    def show(self): print(self.string())

    def string(self): return 'Name: {} Elem: {} X: {} Y: {} Z:{}'.format(self.name, self.elem, self.x, self.y, self.z)

    def simple(self): return '({}, {}, {})'.format(self.x, self.y, self.z)

    def from_origin(self): self.from_center = sqrt((self.x ** 2) + (self.y ** 2) + (self.z ** 2))

class Equation:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    def string(self):
        return '{:6.3f}x + {:6.3f}y + {:6.3f}z + {:6.3f} = 0'.format(self.a, self.b, self.c, self.d)

    def show(self):
        print('{:6.3f}x + {:6.3f}y + {:6.3f}z + {:6.3f} = 0'.format(self.a, self.b, self.c, self.d))

    def func_form(self):
        z = -self.c
        return 'z = {:6.3f}x + {:6.3f}y + {:6.3f}'.format(self.a / z, self.b / z, self.d / z)

    def func_form_tup(self):
        z = -self.c
        return self.a / z, self.b / z, self.d / z

class Struct:
    def __init__(self, pdb_id=None, group=None, title=None, ec_nums=None, records=None, org_name=None, org_sci=None,
                 org_taxid=None, ex_sys=None, ex_sys_taxid=None, organ=None, file_type=None, ec_str=None,
                 atoms=None, planes=None):
        self.pdb_id = pdb_id if pdb_id is not None else ''
        self.group = group if group is not None else ''
        self.title = title if title is not None else ''
        self.ec_nums = ec_nums if ec_nums is not None else []
        self.records = records if records is not None else []
        self.org_name = org_name if org_name is not None else 'No common name'
        self.org_sci = org_sci if org_sci is not None else 'No scientific name'
        self.org_taxid = org_taxid if org_taxid is not None else 'No taxonomy ID'
        self.ex_sys = ex_sys if ex_sys is not None else 'No expression system'
        self.ex_sys_taxid = ex_sys_taxid if ex_sys_taxid is not None else 'No expression system taxonomy ID'
        self.organ = organ if organ is not None else 'No organ'
        self.file_type = file_type if file_type is not None else ''
        self.ec_str = ec_str if ec_str is not None else 'No EC numbers'
        self.atoms = atoms if atoms is not None else []
        self.planes = planes if planes is not None else []

    def get_record(self, atom_name, lig_code, chain):
        for record in self.records:
            if record.atom == atom_name and record.lig_code == lig_code and record.chain == chain: return record
        print('!!!Nothing for atom {} in chain {} of ligand {}'.format(atom_name, chain, lig_code))
        return Record()

    def new_plane(self, atom_names, lig_code, chain):
        atoms = []
        for name in atom_names:
            tmp_rec = self.get_record(name, lig_code, chain)
            if tmp_rec is None: return
            tmp_atm = atom_from_record(tmp_rec)
            tmp_atm.from_origin()
            self.atoms.append(tmp_atm)
            atoms.append(tmp_atm)
        plane = Plane(atoms)
        plane.set_eqn()
        plane.ligand = lig_code
        plane.chain = chain
        self.planes.append(plane)

class Record:
    def __init__(self, pdb_id=DS, label=DS, sn=DS, atom=DS, alt_loc=DS, lig_code=DS, chain=DS, lig_seq=DS, icode=DS,
                 x=DF, y=DF, z=DF, occp=DS, temp=DS, elem=DS, charge=DS):
        self.alt_loc = alt_loc
        self.atom = atom
        self.chain = chain
        self.charge = charge
        self.elem = elem
        self.icode = icode
        self.label = label
        self.lig_code = lig_code
        self.lig_seq = lig_seq
        self.occp = occp
        self.pdb_id = pdb_id
        self.sn = sn
        self.temp = temp
        self.x = x
        self.y = y
        self.z = z

    def info_str(self):
        out = '> {}, {} | Atom: {} | S#: {} | Elem: {} | Lig: {} | LigSeq: {} | X, Y, Z: ({}, {}, {}) | ' \
              'Chain: {}\n'.format(self.pdb_id, self.label, self.atom, self.sn, self.elem, self.lig_code,
                                   self.lig_seq, self.x, self.y, self.z, self.chain)
        return out

class Plane:
    def __init__(self, atoms: [Atom], eqn=None, ligand=DS, chain=DS):
        self.atoms = atoms
        self.eqn = eqn if eqn is not None else Equation(DF, DF, DF, DF)
        self.ligand = ligand
        self.chain = chain

    def string(self):
        atoms = ''
        for atom in self.atoms: atoms += atom.simple() + ' '
        return '{}-{} {}'.format(self.ligand, self.chain, atoms)

    def set_eqn(self):
        if (len(self.atoms)) > 3: self.find_best_fit()
        else: self.eqn = find_plane_eqn(self.atoms[0], self.atoms[1], self.atoms[2])

    def find_best_fit(self):
        x = list((atom.x for atom in self.atoms))
        y = list((atom.y for atom in self.atoms))
        z = list((atom.z for atom in self.atoms))
        names = list((o.name for o in self.atoms))
        modded = [",".join(map(str, comb)) for comb in combinations(names, 3)]
        eqns, modded_list = [], []
        for m in modded: modded_list.append(m.split(','))
        for m in modded_list: eqns.append(find_plane_eqn(self.get_atom(m[0]), self.get_atom(m[1]), self.get_atom(m[2])))
        a, b, c, d = DF, DF, DF, DF
        count = 0

        for eqn in eqns:
            count += 1
            a, b, c, d = a + eqn.a, b + eqn.b, c + eqn.c, d + eqn.d
            if count > 1: a, b, c, d = a / 2, b / 2, c / 2, d / 2

        xyz, avgs = np.array([x, y, z]), np.array([a, b, c, d])
        least_sq_res = leastsq(find_residuals, avgs, args=(None, xyz))[0]
        avg_eqn = Equation(float(least_sq_res[0]), float(least_sq_res[1]), float(least_sq_res[2]),
                           float(least_sq_res[3]))
        self.eqn = avg_eqn

    def get_atom(self, name):
        for atom in self.atoms:
            if name == atom.name: return atom
        print('!!!Could not find atom with name {}, returning blank atom'.format(name))
        return Atom()


def atom_from_record(rec: Record):
    return Atom(rec.atom, rec.elem, rec.x, rec.y, rec.z, ligand=rec.lig_code, chain=rec.chain)

def new_record(line, name):
    """Make new record object from PDB file. Will only create coordinate properties for ATOM/HETATM records."""
    if len(line) < 60: tmp_rec = Record()
    else:
        # ranges taken from the PDB documentation, column values are found in record_formats.txt
        # beginning of range listed in file is one less due to list indices starting at 0. the end value is the same
        # since taking a subsection of a string works as such: substring = string[start:end-1]
        label = skin(line[0:6])
        if K_ATM in label or K_HAT in label:
            x = get_coord(skin(line[30:38]).replace(':', ''))
            y = get_coord(skin(line[38:46]).replace(':', ''))
            z = get_coord(skin(line[46:54]).replace(':', ''))
        else: x, y, z = None, None, None
        tmp_rec = Record(pdb_id=name, label=label, sn=skin(line[6:11]), atom=skin(line[12:16]),
                         alt_loc=line[16:17], lig_code=skin(line[17:20]), chain=skin(line[21:22]),
                         lig_seq=skin(line[22:26]), icode=line[26:27], x=x, y=y, z=z, occp=skin(line[54:60]),
                         temp=skin(line[60:66]), elem=skin(line[76:78]), charge=skin(line[78:80]))
    return tmp_rec

def new_struct(lines, pdb_id):
    struct = Struct()
    struct.pdb_id = pdb_id
    for line in lines:
        if K_CMP in line and K_EC in line:
            ec_nums = re.findall(RE_EC, line)
            for num in ec_nums: struct.ec_nums.append('EC:' + num)
            struct.ec_str = str(struct.ec_nums).replace('[', '').replace(']', '').replace('\'', '')
        if K_TTL in line: struct.title += re.sub(RE_XTRA_SP, ' ', line[10:80])
        if K_HEAD in line and K_REV not in line:
            struct.group += ''.join(re.findall(RE_WORDS, line[10:50]))

        if K_SRC in line:
            try: chunk = re.findall(RE_REC_VAL, line)[0]
            except IndexError: continue

            if K_ORG_CMN in line: struct.org_name = chunk.strip()
            if K_ORG_SCI in line: struct.org_sci = chunk.strip()
            if K_ORG_TAX in line: struct.org_taxid = 'Tax ID: ' + chunk.strip()
            if K_ORG in line: struct.organ = chunk.strip()
            if K_EX_SYS in line: struct.ex_sys = chunk.strip()
            if K_EX_TAX in line: struct.ex_sys_taxid = chunk.strip()

    return struct

def get_coord(coord):
    """Parses string for float coordinate. Blanks, Nones & errors result in return value of 0.0."""
    if coord is None or coord == '': return 0.0
    else:
        tmp_coord = ''
        for c in coord:
            if c == '.' or c == '-' or c.isdigit(): tmp_coord += c
        if tmp_coord == '': return 0.0
        else: return float(tmp_coord)

def get_vector(a: Atom, b: Atom):
    """ Returns (Bx-Ax, By-Ay, Bz-Az) """
    return b.x - a.x, b.y - a.y, b.z - a.z

def find_plane_eqn(i, j, k):
    """Finds components of plane formula ax + by + cz + d = 0 & creates equation object with calculated values."""
    u, v = get_vector(i, j), get_vector(i, k)
    a = (u[1] * v[2]) - (v[1] * u[2])  # a = (By-Ay)(Cz-Az) - (Cy-Ay)(Bz-Az)
    b = (u[2] * v[0]) - (v[2] * u[0])  # b = (Bz-Az)(Cx-Ax) - (Cz-Az)(Bx-Ax)
    c = (u[0] * v[1]) - (v[0] * u[1])  # c = (Bx-Ax)(Cy-Ay) - (Cx-Ax)(By-Ay)
    d = -((a * i.x) + (b * i.y) + (c * i.z))  # d = -(aAx + bAy + cAz)
    eqn = Equation(a, b, c, d)  # make equation object
    return eqn

def f_min(val, p):
    plane = p[0:3]
    dist = (plane * val.T).sum(axis=1) + p[3]
    return dist / np.linalg.norm(plane)

def find_residuals(params, signal, val): return f_min(val, params)

def show_atom_distance(atm1, atm2):
    """Print distance between two atoms with float formatting."""
    dist = atom_distance(atm1, atm2)
    print('distance between {} ({:6.3f}, {:6.3f}, {:6.3f}) & {} ({:6.3f}, {:6.3f}, {:6.3f}) is '
          '~{:6.4f} angstroms'.format(atm1.atom, atm1.x, atm1.y, atm1.z, atm2.atom, atm2.x, atm2.y, atm2.z, dist))

def atom_distance(atm1, atm2):
    """Find distance between two atoms"""
    return sqrt(((atm2.x - atm1.x) ** 2) + ((atm2.y - atm1.y) ** 2) + ((atm2.z - atm1.z) ** 2))

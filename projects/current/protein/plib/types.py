from itertools import combinations
from math import sqrt
from plib.prconstants import *
from plib.prconstants import *
import inspect
import re
import sys

try:
    from scipy.optimize import leastsq
    import numpy as np
except ImportError:
    print('Program needs scipy and numpy to work. One or both are missing. In the terminal, try:\n'
          '`pip3 install scipy` or `pip install scipy`\n'
          '`pip3 install numpy` or `pip install numpy`')
    exit(1)
try: from util import skin
except ImportError:
    # allows for imports from directories at the same level
    sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))
    from sharedlib.util import *

# default values
DS = ' '  # default string
DF = 0.0  # default float
DI = 0  # default int
DL = []  # default list
NS = 'NONE'

class Equation:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def func_form(self): return 'z = {:6.3f}x + {:6.3f}y + {:6.3f}'.format(self.a / -self.c, self.b / -self.c,
                                                                           self.d / -self.c)

    def func_form_tup(self): return self.a / -self.c, self.b / -self.c, self.d / -self.c

    def show(self): print('{:6.3f}x + {:6.3f}y + {:6.3f}z + {:6.3f} = 0'.format(self.a, self.b, self.c, self.d))

    def string(self): return '{:6.3f}x + {:6.3f}y + {:6.3f}z + {:6.3f} = 0'.format(self.a, self.b, self.c, self.d)

class AtomRec:
    def __init__(self, pdb_id=DS, label=DS, sn=DS, atom=DS, alt_loc=DS, lig_code=DS, chain=DS, seq=DS, icode=DS, x=DF,
                 y=DF, z=DF, occp=DS,
                 temp=DS, elem=DS, charge=DS):
        self.alt_loc = alt_loc
        self.atom = atom
        self.chain = chain
        self.charge = charge
        self.elem = elem
        self.icode = icode
        self.label = label
        self.lig_code = lig_code
        self.occp = occp
        self.pdb_id = pdb_id
        self.seq = seq
        self.sn = sn
        self.temp = temp
        self.x = x
        self.y = y
        self.z = z

    def info_str(self):
        out = '> {}, {} | ATOM {} | SN {} | ELEM {} | LIG {} | SEQ {} | X, Y, Z: ({}, {}, {}) | ' \
              'CHAIN {}\n'.format(self.pdb_id, self.label, self.atom, self.sn, self.elem, self.lig_code,
                                  self.seq, self.x, self.y, self.z, self.chain)
        return out

    def simple(self): return '{}=({},{},{})'.format(self.atom, self.x, self.y, self.z)

    def show(self): print(
        '{}: {} {} SEQ_{} CHAIN_{}'.format(self.pdb_id, self.atom, self.lig_code, self.seq, self.chain))

class Ring:
    def __init__(self, atoms: [AtomRec], eqn=None, ligand=DS, chain=DS, seq=DS, ring_type=DS):
        self.atoms = atoms
        self.chain = chain
        self.eqn = eqn if eqn is not None else Equation(DF, DF, DF, DF)
        self.ligand = ligand
        self.seq = seq
        self.ring_type = ring_type

    def string(self):
        atoms = ''
        for atom in self.atoms: atoms += atom.simple() + ' '
        return '{}: SEQ_{} RING_{} CHAIN_{} {}'.format(self.ligand, self.seq, self.ring_type, self.chain, atoms)

    def set_eqn(self):
        if (len(self.atoms)) > 3: self.find_best_fit()
        else:
            if len(self.atoms) == 3: self.eqn = find_plane_eqn(self.atoms[0], self.atoms[1], self.atoms[2])
            else: return Equation(0.0, 0.0, 0.0, 0.0)

    def find_best_fit(self):
        x = list((atom.x for atom in self.atoms))
        y = list((atom.y for atom in self.atoms))
        z = list((atom.z for atom in self.atoms))
        names = list((o.atom for o in self.atoms))
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
            if name == atom.atom: return atom
        print('!!!Could not find atom with name {}, returning blank atom'.format(name))
        return Atom()

class ProtStruct:
    def __init__(self, pdb_id=NS, group=NS, title=NS, ec_nums=DL, records: [AtomRec] = DL, org_name=NS, org_sci=NS,
                 org_taxid=NS, ex_sys=NS,
                 ex_sys_taxid=NS, organ=NS, ft=NS, ec_str=NS, atoms: [AtomRec] = DL, rings: [Ring] = DL, ligands=DL):
        self.ec_nums = ec_nums
        self.ec_str = ec_str
        self.ex_sys = ex_sys
        self.ex_sys_taxid = ex_sys_taxid
        self.ft = ft
        self.group = group
        self.ligands = ligands
        self.org_name = org_name
        self.org_sci = org_sci
        self.org_taxid = org_taxid
        self.organ = organ
        self.pdb_id = pdb_id
        self.rings = rings
        self.records = records
        self.title = title

    def entry_str(self, num_rec):
        return self.pdb_id + ' || ' + self.group + ' || ' + self.org_name + ' || ' + self.org_taxid + ' || ' + self.ex_sys + ' || ' + \
               (str(self.ec_nums)) + ' || ' + str(num_rec) + ' records'

    def get_record(self, atom_name, lig_code, chain, seq):
        for record in self.records:
            try:
                if record.atom == atom_name and record.lig_code == lig_code and record.chain == chain and seq == record.seq: return record
            except AttributeError: return None
            # print('record was type {}, function called from {}. object was ring {}'.format(type(record), inspect.stack()[1][3],
            #                                                                                record.string()))
            # except AttributeError: pass
            # print('!!!Nothing for atom {} in chain {} of ligand {}'.format(atom_name, chain, lig_code))

    def new_ring(self, atom_names, lig_code, chain, seq, ring_name):
        atoms = []
        for name in atom_names:
            atom_rec = self.get_record(name, lig_code, chain, seq)
            if atom_rec is None:
                # print('did not find atom record, function called from {}'.format(inspect.stack()[1][3]))
                return
            atoms.append(atom_rec)
        ring = Ring(atoms=atoms, ligand=lig_code, chain=chain, seq=seq, ring_type=ring_name)
        ring.set_eqn()
        self.rings.append(ring)

def new_record(line, name):
    """Make new record object from PDB file. Will only create coordinate properties for ATOM/HETATM records."""
    if len(line) < 60: tmp_rec = AtomRec(pdb_id=name)
    else:
        # ranges taken from the PDB documentation, column values are found in record_formats.txt beginning of range listed in file is one
        # less due to list indices starting at 0. the end value is the same since taking a subsection of a string works as such:
        # substring = string[start:end-1]
        label = skin(line[0:6])
        if K_ATM in label or K_HAT in label:
            coord_line = re.findall(r'(-?\d+\.*\d*\s*)(-?\d+\.*\d*\s*)(-?\d+\.*\d*\s)', line[30:54])
            if not len(coord_line): x, y, z = None, None, None
            else:
                match = coord_line[0]
                x = get_coord(skin(match[0]).replace(':', ''), name)
                y = get_coord(skin(match[1]).replace(':', ''), name)
                z = get_coord(skin(match[2]).replace(':', ''), name)
        else: x, y, z = None, None, None
        tmp_rec = AtomRec(pdb_id=name, label=label, sn=skin(line[6:11]), atom=skin(line[12:16]),
                          alt_loc=line[16:17], lig_code=skin(line[17:20]), chain=skin(line[21:22]),
                          seq=skin(line[22:26]), icode=line[26:27], x=x, y=y, z=z, occp=skin(line[54:60]),
                          temp=skin(line[60:66]), elem=skin(line[76:78]), charge=skin(line[78:80]))
    return tmp_rec

def new_struct(lines, pdb_id):
    struct = ProtStruct()
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

def get_coord(coord, pdb_id):
    """Parses string for float coordinate. Blanks, Nones & errors result in return value of 0.0."""
    if coord is None or coord == '': return 0.0
    else:
        tmp_coord = ''
        for c in coord:
            if c == '.' or coord.index(c) or c.isdigit(): tmp_coord += c
        try:
            if tmp_coord == '': return 0.0
            else: return float(tmp_coord)
        except ValueError as e:
            print('error with file {}: {}. input val = {}'.format(pdb_id, e, coord))
            return 0.0

def get_vector(a: AtomRec, b: AtomRec):
    """ Returns (Bx-Ax, By-Ay, Bz-Az) """
    return b.x - a.x, b.y - a.y, b.z - a.z

def find_plane_eqn(pa, pb, pc):
    """Finds components of plane formula ax + by + cz + d = 0 & creates equation object with calculated values."""
    ab, ac = get_vector(pa, pb), get_vector(pa, pc)
    a = (ab[1] * ac[2]) - (ac[1] * ab[2])  # a = (By-Ay)(Cz-Az) - (Cy-Ay)(Bz-Az)
    b = (ab[2] * ac[0]) - (ac[2] * ab[0])  # b = (Bz-Az)(Cx-Ax) - (Cz-Az)(Bx-Ax)
    c = (ab[0] * ac[1]) - (ac[0] * ab[1])  # c = (Bx-Ax)(Cy-Ay) - (Cx-Ax)(By-Ay)
    d = -((a * pa.x) + (b * pa.y) + (c * pa.z))  # d = -(aAx + bAy + cAz)
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

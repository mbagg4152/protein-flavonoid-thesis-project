import re
from plib.strings_consts import *
import sys

sys.path.append(os.getcwd().replace(os.sep + 'protein', ''))  # allows for imports from directories at the same level
from lib.jsondata import *
from lib.util import *

class Record:
    def __init__(self, pdb_id=None, label=None, sn=None, atom=None, alt_loc=None, lig_code=None, chain_id=None,
                 lig_seq=None, ins_code=None, x=None, y=None, z=None, occupy=None, temp=None, elem=None, charge=None):
        self.pdb_id = pdb_id if pdb_id is not None else ' '
        self.label = label if label is not None else ' '
        self.sn = sn if sn is not None else ' '
        self.atom = atom if atom is not None else ' '
        self.alt_loc = alt_loc if alt_loc is not None else ' '
        self.lig_code = lig_code if lig_code is not None else ' '
        self.chain_id = chain_id if chain_id is not None else ' '
        self.lig_seq = lig_seq if lig_seq is not None else ' '
        self.ins_code = ins_code if ins_code is not None else ' '
        self.x = x if x is not None else 0.0
        self.y = y if y is not None else 0.0
        self.z = z if z is not None else 0.0
        self.occupy = occupy if occupy is not None else ' '
        self.temp = temp if temp is not None else ' '
        self.elem = elem if elem is not None else ' '
        self.charge = charge if charge is not None else ' '

    def info_str(self):
        out = '> {}, {} | Atom: {} | Serial #: {} | Elem: {} | Lig: {} | LigSeq: {} | X, Y, Z: ({}, {}, {}) | ' \
              'Chain: {}\n'.format(self.pdb_id, self.label, self.atom, self.sn, self.elem, self.lig_code,
                                   self.lig_seq, self.x, self.y, self.z, self.chain_id)
        return out


class Atom:
    def __init__(self, name=None, elem=None, x=None, y=None, z=None):
        self.name = name if name is not None else ' '
        self.elem = elem if elem is not None else ' '
        self.x = x if x is not None else 0.0
        self.y = y if y is not None else 0.0
        self.z = z if z is not None else 0.0

    def show(self):
        print(self.string())

    def string(self):
        return 'Name: {} Elem: {} X: {} Y: {} Z:{}'.format(self.name, self.elem, self.x, self.y, self.z)

    def simple(self):
        return '[{} ({}, {}, {})]'.format(self.name, self.x, self.y, self.z)

class Plane:
    def __init__(self, atoms: [Atom]):
        self.atoms = atoms

    def string(self):
        atoms = ''
        for atom in self.atoms: atoms += atom.simple() + ' '
        return 'Plane: {}'.format(atoms)


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

    def get_record(self, atom_name):
        for record in self.records:
            if record.atom == atom_name: return record
        print('!!!WARNING!!!No record found containing atom name {}! Returning empty record'.format(atom_name))
        return Record()

    def new_plane(self, atom_names):
        atoms = []
        for name in atom_names:
            tmp_rec = self.get_record(name)
            tmp_atm = new_atom(tmp_rec)
            self.atoms.append(tmp_atm)
            atoms.append(tmp_atm)
        plane = Plane(atoms)
        self.planes.append(plane)

def new_atom(rec: Record):
    atom = Atom()
    atom.name = rec.atom
    atom.elem = rec.elem
    atom.x = rec.x
    atom.y = rec.y
    atom.z = rec.z
    return atom

def new_record(line, name):
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
        else:
            x = None
            y = None
            z = None
        tmp_rec = Record(pdb_id=name, label=label, sn=skin(line[6:11]), atom=skin(line[12:16]),
                         alt_loc=line[16:17], lig_code=skin(line[17:20]), chain_id=skin(line[21:22]),
                         lig_seq=skin(line[22:26]), ins_code=line[26:27], x=x, y=y, z=z, occupy=skin(line[54:60]),
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
    if coord is None or coord == '': return 0.0
    else:
        tmp_coord = ''
        # print('coord val: ' + coord)
        for c in coord:
            if c == '.' or c == '-' or c.isdigit(): tmp_coord += c
        if tmp_coord == '': return 0.0
        else: return float(tmp_coord)

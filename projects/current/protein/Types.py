import re


class Record:
    def __init__(self, pdb_id, label, serial, atom_name, alt_loc, ligand_code, chain_id, ligand_seq, ins_code, pos_x,
                 pos_y, pos_z, occupy, temp, elem, charge):
        self.pdb_id = pdb_id if pdb_id is not None else ''
        self.label = label if label is not None else ''
        self.serial = serial if serial is not None else ''
        self.atom_name = atom_name if atom_name is not None else ''
        self.alt_loc = alt_loc if alt_loc is not None else ''
        self.ligand_code = ligand_code if ligand_code is not None else ''
        self.chain_id = chain_id if chain_id is not None else ''
        self.ligand_seq = ligand_seq if ligand_seq is not None else ''
        self.ins_code = ins_code if ins_code is not None else ''
        self.pos_x = pos_x if pos_x is not None else ''
        self.pos_y = pos_y if pos_y is not None else ''
        self.pos_z = pos_z if pos_z is not None else ''
        self.occupy = occupy if occupy is not None else ''
        self.temp = temp if temp is not None else ''
        self.elem = elem if elem is not None else ''
        self.charge = charge if charge is not None else ''

    def important_str(self):
        out = '> ' + self.pdb_id + ', ' + self.label + ' | AtomName: ' + self.atom_name + ' | Serial: ' + self.serial \
              + '| Elem: ' + self.elem + ' | Ligand: ' + self.ligand_code + ' | LigandSeq: ' + self.ligand_seq + \
              ' | X,Y,Z: (' + self.pos_x + ', ' + self.pos_y + ', ' + self.pos_z + ') | ChainID: ' + self.chain_id + \
              '\n'
        return out


class Entry:
    def __init__(self, pdb_id=None, header=None, title=None, ec_nums=None, records=None):
        self.pdb_id = pdb_id if pdb_id is not None else ''
        self.header = header if header is not None else ''
        self.title = title if title is not None else ''
        self.ec_nums = ec_nums if ec_nums is not None else []
        self.records = records if records is not None else []


def new_record(line, name):
    if len(line) < 60: return None
    else:
        # ranges taken from the PDB documentation, column values are found in record_formats.txt
        # beginning of range listed in file is one less due to list indices starting at 0. the end value is the same
        # since taking a subsection of a string works as such: substring = string[start:end-1]
        tmp_record = Record(pdb_id=name, label=line[0:6], serial=line[6:11], atom_name=line[12:16], alt_loc=line[16:17],
                            ligand_code=line[17:20], chain_id=line[21:22], ligand_seq=line[22:26], ins_code=line[26:27],
                            pos_x=line[30:38], pos_y=line[38:45], pos_z=line[46:54], occupy=line[54:60],
                            temp=line[60:66], elem=line[76:78], charge=line[78:80])
        return tmp_record


def new_entry(lines):
    tmp_entry = Entry()
    excess_space = ' {2,}'
    for line in lines:
        if 'COMPND' in line and 'EC:' in line:
            ec_pattern = '([0-9]\.{1}[^ ,;]*)'
            tmp_entry.ec_nums = re.findall(ec_pattern, line)
            pass
        if 'TITLE' in line: tmp_entry.title += re.sub(excess_space, ' ', line[10:80])
        if 'HEADER' in line:
            tmp_entry.header += re.sub(excess_space, ' ', line[10:50])
            tmp_entry.pdb_id += re.sub(excess_space, ' ', line[62:66])

    return tmp_entry

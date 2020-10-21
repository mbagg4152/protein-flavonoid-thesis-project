import re
from StringsAndConsts import *


class Record:
    def __init__(self, pdb_id=None, label=None, serial=None, atom_name=None, alt_loc=None, ligand_code=None,
                 chain_id=None, ligand_seq=None, ins_code=None, pos_x=None, pos_y=None, pos_z=None, occupy=None,
                 temp=None, elem=None, charge=None, file_type=None):
        self.pdb_id = pdb_id if pdb_id is not None else ' '
        self.label = label if label is not None else ' '
        self.serial = serial if serial is not None else ' '
        self.atom_name = atom_name if atom_name is not None else ' '
        self.alt_loc = alt_loc if alt_loc is not None else ' '
        self.ligand_code = ligand_code if ligand_code is not None else ' '
        self.chain_id = chain_id if chain_id is not None else ' '
        self.ligand_seq = ligand_seq if ligand_seq is not None else ' '
        self.ins_code = ins_code if ins_code is not None else ' '
        self.pos_x = pos_x if pos_x is not None else ' '
        self.pos_y = pos_y if pos_y is not None else ' '
        self.pos_z = pos_z if pos_z is not None else ' '
        self.occupy = occupy if occupy is not None else ' '
        self.temp = temp if temp is not None else ' '
        self.elem = elem if elem is not None else ' '
        self.charge = charge if charge is not None else ' '
        self.file_type = file_type if file_type is not None else ' '

    def important_str(self):
        out = '> ' + self.pdb_id + ', ' + self.label + ' | AtomName: ' + self.atom_name + ' | Serial: ' + self.serial \
              + '  | Elem: ' + self.elem + ' | Ligand: ' + self.ligand_code + ' | LigandSeq: ' + self.ligand_seq + \
              ' | X,Y,Z: (' + self.pos_x + ', ' + self.pos_y + ', ' + self.pos_z + ') | ChainID: ' + self.chain_id + \
              '\n'
        return out


class Entry:
    def __init__(self, pdb_id=None, group=None, title=None, ec_nums=None, records=None, org_name=None, org_sci=None,
                 org_taxid=None, ex_sys=None, ex_sys_taxid=None, organ=None, file_type=None):
        self.pdb_id = pdb_id if pdb_id is not None else ' '
        self.group = group if group is not None else ' '
        self.title = title if title is not None else ' '
        self.ec_nums = ec_nums if ec_nums is not None else []
        self.records = records if records is not None else []
        self.org_name = org_name if org_name is not None else ' '
        self.org_sci = org_sci if org_sci is not None else ' '
        self.org_taxid = org_taxid if org_taxid is not None else ' '
        self.ex_sys = ex_sys if ex_sys is not None else ' '
        self.ex_sys_taxid = ex_sys_taxid if ex_sys_taxid is not None else ' '
        self.organ = organ if organ is not None else ' '
        self.file_type = file_type if file_type is not None else ' '


def new_record(line, name, file_type):
    if len(line) < 60: tmp_record = Record()
    else:
        # ranges taken from the PDB documentation, column values are found in record_formats.txt
        # beginning of range listed in file is one less due to list indices starting at 0. the end value is the same
        # since taking a subsection of a string works as such: substring = string[start:end-1]
        if file_type == 'PDB':
            tmp_record = Record(pdb_id=name, label=line[0:6], serial=line[6:11], atom_name=line[12:16],
                                alt_loc=line[16:17], ligand_code=line[17:20], chain_id=line[21:22],
                                ligand_seq=line[22:26], ins_code=line[26:27], pos_x=line[30:38], pos_y=line[38:45],
                                pos_z=line[46:54], occupy=line[54:60], temp=line[60:66], elem=line[76:78],
                                charge=line[78:80], file_type=file_type)
        else:
            tmp_record = Record()
    return tmp_record


def new_entry(lines):
    tmp_entry = Entry()

    for line in lines:
        if K_CMP in line and K_EC in line:
            tmp_entry.ec_nums = re.findall(RE_EC, line)
            pass
        if K_TTL in line: tmp_entry.title += re.sub(RE_XTRA_SP, ' ', line[10:80])
        if K_HEAD in line:
            tmp_entry.group += re.sub(RE_XTRA_SP, ' ', line[10:50])
            tmp_entry.pdb_id += re.sub(RE_XTRA_SP, ' ', line[62:66])

        if K_SRC in line:
            try: chunk = re.findall(RE_REC_VAL, line)[0]
            except IndexError: chunk = ''

            if K_ORG_CMN in line: tmp_entry.org_name = chunk
            if K_ORG_SCI in line: tmp_entry.org_sci = chunk
            if K_ORG_TAX in line: tmp_entry.org_taxid = chunk
            if K_ORG in line: tmp_entry.organ = chunk
            if K_EX_SYS in line: tmp_entry.ex_sys = chunk
            if K_EX_TAX in line: tmp_entry.ex_sys_taxid = chunk
    return tmp_entry

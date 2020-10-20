class Record:
    def __init__(self, pdb_id, label, serial, atom_name, alt_loc, ligand_code, chain_id, ligand_seq, ins_code, pos_x,
                 pos_y, pos_z, occupy, temp, elem, charge):
        self.pdb_id = pdb_id
        self.label = label
        self.serial = serial
        self.atom_name = atom_name
        self.alt_loc = alt_loc
        self.ligand_code = ligand_code
        self.chain_id = chain_id
        self.ligand_seq = ligand_seq
        self.ins_code = ins_code
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.pos_z = pos_z
        self.occupy = occupy
        self.temp = temp
        self.elem = elem
        self.charge = charge

    def important_str(self):
        out = '> ' + self.pdb_id + ', ' + self.label + ' | AtomName: ' + self.atom_name + ' | Serial: ' + self.serial \
              + '| Elem: ' + self.elem + ' | Ligand: ' + self.ligand_code + ' | LigandSeq: ' + self.ligand_seq + \
              ' | X,Y,Z: (' + self.pos_x + ', ' + self.pos_y + ', ' + self.pos_z + ') | ChainID: ' + self.chain_id + \
              '\n'
        return out


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

REC = (0, 5)
S_NUM = (6, 10)
A_NAME = (12, 15)
ALT = (16, 16)
RES_NAME = (17, 19)
CHAIN_ID = (21, 21)
RES_SEQ = (22, 25)
I_CODE = (26, 26)
ANG_X = (30, 37)
ANG_Y = (38, 45)
ANG_Z = (46, 53)
OCC = (54, 59)
TEMP = (60, 65)
ELEM = (76, 78)
CHARGE = (78, 79)


class Record:
    def __init__(self, label, serial, atom, alt, lig, chain, lig_seq, ic, x, y, z, occ, tf, elem, charge):
        self.label = label
        self.serial = serial
        self.atom = atom
        self.alt_loc = alt
        self.lig_code = lig
        self.chain_id = chain
        self.lig_seq = lig_seq
        self.icode = ic
        self.x = x
        self.y = y
        self.z = z
        self.occupy = occ
        self.temp = tf
        self.elem = elem
        self.charge = charge

    def full_str(self):
        out = 'Record: ' + self.label + ' AtomSN: ' + self.serial + ' AtomName: ' + self.atom + \
              ' AltLoc: ' + self.alt_loc + ' ResName: ' + self.lig_code + ' ChainID: ' + \
              self.chain_id + ' ResSeq: ' + self.lig_seq + ' InsCode: ' + self.icode + ' X: ' + \
              self.x + ' Y: ' + self.y + ' Z: ' + self.z + ' Occupancy: ' + self.occupy + ' temp factor' + self.temp + \
              ' Elem: ' + self.elem + ' Charge: ' + self.charge + '\n'
        return out

    def relevant_str(self):
        out = self.label + ' Elem: ' + self.elem + ' AtomSerial: ' + self.serial + ' AtomName: ' + self.atom \
              + ' Ligand: ' + self.lig_code + ' ChainID: ' + self.chain_id + ' LigandSeq: ' + self.lig_seq + \
              ' IC: ' + self.icode + ' X: ' + self.x + ' Y: ' + self.y + ' Z: ' + self.z + '\n'
        return out
def new_record(line):
    if len(line) < 60: return None
    else:
        rec = line[REC[0]:REC[1] + 1]
        serial = line[S_NUM[0]:S_NUM[1] + 1]
        atom = line[A_NAME[0]:A_NAME[1] + 1]
        alt = line[ALT[0]:ALT[1] + 1]
        res_name = line[RES_NAME[0]:RES_NAME[1] + 1]
        chain = line[CHAIN_ID[0]:CHAIN_ID[1] + 1]
        res_seq = line[RES_SEQ[0]:RES_SEQ[1] + 1]
        icode = line[I_CODE[0]:I_CODE[1] + 1]
        x = line[ANG_X[0]:ANG_X[1] + 1]
        y = line[ANG_Y[0]:ANG_Y[1] + 1]
        z = line[ANG_Z[0]:ANG_Z[1] + 1]
        occ = line[OCC[0]:OCC[1] + 1]
        tf = line[TEMP[0]:TEMP[1] + 1]
        elem = line[ELEM[0]:ELEM[1] + 1]
        charge = line[CHARGE[0]:CHARGE[1] + 1]
        tmp_record = Record(rec, serial, atom, alt, res_name, chain, res_seq, icode, x, y, z, occ, tf, elem, charge)
        return tmp_record

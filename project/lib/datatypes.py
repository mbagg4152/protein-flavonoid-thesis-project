from lib.pathstrings import *
from lib.miscvals import *

class ChemData:
    def __init__(self, label: str, species: [str], file_name: str):
        self.species = species
        self.label = label
        self.file_name = file_name

class Species:
    def __init__(self, name: str, count: int, flavonoids: [str]):
        self.name = name
        self.flavonoids = flavonoids
        self.count = count

    def species_string(self):
        flav_str = ''
        for i, f in enumerate(self.flavonoids):
            flav_str += f
            if i != len(self.flavonoids) - 1:
                flav_str += ', '
        line = '' + str(self.count) + ', ' + self.name + ', ' + flav_str
        return line


cd_api = ChemData(label=AGI, species=[], file_name=FN_AGI)
cd_bun = ChemData(label=BUN, species=[], file_name=FN_BUN)
cd_ec = ChemData(label=EC, species=[], file_name=FN_EC)
cd_egt = ChemData(label=EGT, species=[], file_name=FN_EGT)
cd_erd = ChemData(label=ERD, species=[], file_name=FN_ERD)
cd_gc = ChemData(label=GC, species=[], file_name=FN_GC)
cd_gen = ChemData(label=GEN, species=[], file_name=FN_GEN)
cd_hwb = ChemData(label=HWB, species=[], file_name=FN_HWB)
cd_kmp = ChemData(label=KMP, species=[], file_name=FN_KMP)
cd_kxn = ChemData(label=KXN, species=[], file_name=FN_KXN)
cd_lu2 = ChemData(label=LU2, species=[], file_name=FN_LU2)
cd_myc = ChemData(label=MYC, species=[], file_name=FN_MYC)
cd_nar = ChemData(label=NAR, species=[], file_name=FN_NAR)
cd_que = ChemData(label=QUE, species=[], file_name=FN_QUER)
cd_hcc = ChemData(label=HCC, species=[], file_name=FN_HCC)
data_lists = [cd_api, cd_bun, cd_kxn, cd_hwb, cd_ec, cd_egt, cd_erd, cd_gc, cd_gen, cd_kmp, cd_lu2, cd_myc, cd_nar,
              cd_que, cd_hcc]

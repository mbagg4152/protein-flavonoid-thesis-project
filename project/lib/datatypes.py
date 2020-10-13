from lib.pathstrings import *
from lib.compoundinfo import *
from lib.miscstrings import *

class ChemData:
    def __init__(self, label: str, species: [str], logic, file_name: str):
        self.logic = logic
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

class ThreadRes:
    def __init__(self, path_id, gene_data):
        self.path_id = path_id
        self.gene_data = gene_data

class Flavonoid:
    def __init__(self):
        self.name = ''
        self.req = ''
        self.plants = []

class Plant:
    def __init__(self):
        self.name = ''
        self.code = ''
        self.path_codes = []
        self.ec_numbers = []
        self.flavonoids = []

out_data = [(AGI, [], agi, fn_apig), (BU, [], bu, fn_bute), (CA, [], ca, fn_cate),
            (CYN, [], cyn, fn_cyan), (EC, [], ecat, fn_ecat), (EGCG, [], egcg, fn_epig),
            (ERD, [], erd, fn_erio), (GC, [], gc, fn_gall), (GEN, [], gen, fn_geni),
            (KMP, [], kmp, fn_kaem), (LU2, [], lu2, fn_lute), (MYC, [], myri, fn_myri),
            (NAR, [], nar, fn_nari), (QUE, [], quer, fn_quer), (HCC, [], bu, fn_isol)]
# out_data = [(AGI, [], agi, fn_apig), (BU, [], bu, fn_bute), (CA, [], ca, fn_cate),
#             (CYN, [], cyn, fn_cyan), (EC, [], ecat, fn_ecat), (EGCG, [], egcg, fn_epig),
#             (ERD, [], erd, fn_erio), (GC, [], gc, fn_gall), (GEN, [], gen, fn_geni),
#             (KMP, [], kmp, fn_kaem), (LU2, [], lu2, fn_lute), (MYC, [], myri, fn_myri),
#             (NAR, [], nar, fn_nari), (QUE, [], quer, fn_quer), (E23, [], ec23, fn_ec23),
#             (E24, [], ec24, fn_ec24), (E25, [], ec25, fn_ec25),
#             (E26, [], ec26, fn_ec26), (HCC, [], bu, fn_isol)]

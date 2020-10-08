from lib.pathstrings import *
from lib.compoundinfo import *

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

out_data = [(APIG, [], has_apig, fn_apig), (BUTE, [], has_bute, fn_bute), (CATE, [], has_cate, fn_cate),
            (CYAN, [], has_cyan, fn_cyan), (ECAT, [], has_ecat, fn_ecat), (EPIG, [], has_epig, fn_epig),
            (ERIO, [], has_erio, fn_erio), (GALL, [], has_gall, fn_gall), (GENI, [], has_geni, fn_geni),
            (KAEM, [], has_kaem, fn_kaem), (LUTE, [], has_lute, fn_lute), (MYRI, [], has_myri, fn_myri),
            (NARI, [], has_nari, fn_nari), (QUER, [], has_quer, fn_quer), (ENW, [], has_w, fn_enz_w),
            (ENX, [], has_x, fn_enz_x), (ENY, [], has_y, fn_enz_y),
            (ENZ, [], has_z, fn_enz_z)]

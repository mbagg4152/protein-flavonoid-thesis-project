from lib.file_path_names import *
from lib.flavonoid_data import *


class ChemData:
    def __init__(self, label: str, species: [str], logic: bool, file_name: str):
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


out_data = [(APIG, [], logic_apig, fn_apig), (BUTE, [], logic_bute, fn_bute), (CATE, [], logic_cate, fn_cate),
            (CYAN, [], logic_cyan, fn_cyan), (ECAT, [], logic_ecat, fn_ecat), (EPIG, [], logic_epig, fn_epig),
            (ERIO, [], logic_erio, fn_erio), (GALL, [], logic_gall, fn_gall), (GENI, [], logic_geni, fn_geni),
            (KAEM, [], logic_kaem, fn_kaem), (LUTE, [], logic_lute, fn_lute), (MYRI, [], logic_myri, fn_myri),
            (NARI, [], logic_nari, fn_nari), (QUER, [], logic_quer, fn_quer), (ENZ_W, [], ENZ_W, fn_enz_w),
            (ENZ_X, [], ENZ_X, fn_enz_x), (ENZ_Y, [], ENZ_Y, fn_enz_3),
            (ENZ_Z, [], ENZ_Z, fn_enz_4)]

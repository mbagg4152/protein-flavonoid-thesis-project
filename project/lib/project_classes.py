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


out_data = [(APIG, [], logic_apig, file_apig), (BUTE, [], logic_bute, file_bute), (CATE, [], logic_cate, file_cate),
            (CYAN, [], logic_cyan, file_cyan), (ECAT, [], logic_ecat, file_ecat), (EPIG, [], logic_epig, file_epig),
            (ERIO, [], logic_erio, file_erio), (GALL, [], logic_gall, file_gall), (GENI, [], logic_geni, file_geni),
            (KAEM, [], logic_kaem, file_kaem), (LUTE, [], logic_lute, file_lute), (MYRI, [], logic_myri, file_myri),
            (NARI, [], logic_nari, file_nari), (QUER, [], logic_quer, file_quer), (logic_ec1, [], logic_ec1, file_ec1),
            (logic_ec2, [], logic_ec2, file_ec2), (logic_ec3, [], logic_ec3, file_ec3),
            (logic_ec4, [], logic_ec4, file_ec4)]

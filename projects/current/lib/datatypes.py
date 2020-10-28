from lib.pathstrings import *
from lib.miscvals import *


class ChemData:
    def __init__(self, label: str, plants: [str], file_name: str):
        self.plants = plants
        self.label = label
        self.file_name = file_name

    def __eq__(self, other):
        return self.plants == other.plants and \
               self.label == other.label and \
               self.file_name == other.file_name

    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False


class EcFastaCollection:
    def __init__(self, ec_num=None, ec_entries=None):
        self.ec_name = ec_num if ec_num is not None else ' '
        self.ec_entries = ec_entries if ec_entries is not None else []

    def __eq__(self, other):
        return self.ec_name == other.ec_name and \
               self.ec_entries == other.ec_entries

    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False


class EcCounts:
    def __init__(self, number=None, count=None):
        self.number = number if number is not None else ' '
        self.count = count if number is not None else 0


class FastaEcEntry:
    def __init__(self, gene=None, dna=None, plant=None):
        self.gene_id = gene if gene is not None else ' '
        self.dna_seq = dna if dna is not None else ' '
        self.plant = plant if dna is not None else ' '

    def __eq__(self, other):
        return self.gene_id == other.gene_id and \
               self.dna_seq == other.dna_seq and \
               self.plant == other.plant

    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False

    def simple(self): return self.dna_seq


class Gene:
    def __init__(self, gene_id=None, plant=None, compound=None, ec_num=None, k_ortho=None, path=None, plant_code=None):
        self.gene_id = gene_id if gene_id is not None else ' '
        self.plant = plant if plant is not None else ' '
        self.plant_code = plant_code if plant_code is not None else ' '
        self.compound = compound if compound is not None else ' '
        self.ec_num = ec_num if ec_num is not None else ' '
        self.k_ortho = k_ortho if k_ortho is not None else ' '
        self.path = path if path is not None else ' '

    def __eq__(self, other):
        return self.gene_id == other.gene_id and \
               self.plant_code == other.plant_code and \
               self.ec_num == other.ec_num

    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False

    def simple(self):
        return self.plant + ', ' + self.gene_id + ', ' + self.compound + ', ' + self.ec_num + ', ' + self.k_ortho

    def no_plant(self):
        return self.gene_id + ', ' + self.compound + ', ' + self.ec_num + ', ' + self.k_ortho


class PathGene:
    def __init__(self, path=None, genes=None):
        self.path = path if path is not None else ' '
        self.genes = genes if genes is not None else []

    def __eq__(self, other):
        return self.path == other.path and \
               self.genes == other.genes

    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False


class Plant:
    def __init__(self, name=None, code=None, genes=None, ec_nums=None, flavonoids=None, ec_counts=None):
        self.name = name if name is not None else ' '
        self.code = code if code is not None else ' '
        self.genes = genes if genes is not None else []
        self.ec_nums = ec_nums if ec_nums is not None else []
        self.flavonoids = flavonoids if flavonoids is not None else []
        self.ec_counts = ec_counts if ec_counts is not None else []

    def __eq__(self, other):
        return self.name == other.name and \
               self.code == other.code and \
               self.genes == other.genes and \
               self.ec_nums == other.ec_nums and \
               self.flavonoids == other.flavonoids

    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False

    def simple(self):
        gstr = ''
        for gene in self.genes: gstr += gene.no_plant() + ' || '
        return self.name + ', ' + self.code + ', ' + gstr + ', ' + str(self.ec_nums) + ', ' + str(self.flavonoids)

    def has_ec_count(self, ec_number):
        for ec in self.ec_counts:
            if ec.number == ec_number: return True
        return False

    def incr_ec_count(self, ec_number):
        for ec in self.ec_counts:
            if ec.number == ec_number:
                ec.count += 1


cd_api = ChemData(label=AGI, plants=[], file_name=FN_AGI)  # apigenin
cd_bun = ChemData(label=BUN, plants=[], file_name=FN_BUN)  # butein
cd_ec = ChemData(label=EC, plants=[], file_name=FN_EC)  # epicatechin
cd_egt = ChemData(label=EGT, plants=[], file_name=FN_EGT)  # epigallocatechin
cd_erd = ChemData(label=ERD, plants=[], file_name=FN_ERD)  # eriodictyol
cd_gc = ChemData(label=GC, plants=[], file_name=FN_GC)  # gallocatechin
cd_gen = ChemData(label=GEN, plants=[], file_name=FN_GEN)  # genistein
cd_hwb = ChemData(label=HWB, plants=[], file_name=FN_HWB)  # cyanidin
cd_kmp = ChemData(label=KMP, plants=[], file_name=FN_KMP)  # kaempferol
cd_kxn = ChemData(label=KXN, plants=[], file_name=FN_KXN)  # catechin
cd_lu2 = ChemData(label=LU2, plants=[], file_name=FN_LU2)  # luteolin
cd_myc = ChemData(label=MYC, plants=[], file_name=FN_MYC)  # myricetin
cd_nar = ChemData(label=NAR, plants=[], file_name=FN_NAR)  # naringenin
cd_que = ChemData(label=QUE, plants=[], file_name=FN_QUER)  # quercetin
cd_hcc = ChemData(label=HCC, plants=[], file_name=FN_HCC)  # isoliquiritigenin
data_lists = [cd_api, cd_bun, cd_kxn, cd_hwb, cd_ec, cd_egt,
              cd_erd, cd_gc, cd_gen, cd_kmp, cd_lu2, cd_myc,
              cd_nar, cd_que, cd_hcc]

from const_vals import *
from const_paths import *


class ChemData:
    """
    This class holds the data for each flavonoid. The objects are initialized with their file name and label and only
    later in the program, their empty list of plants will be filled.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.label: string that contains the flavonoids name
    self.plants: list of plants predicted to produce the flavonoid
    self.file_name: string that holds the flavonoids output file name
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    """

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


class Plant:
    """
    This object holds information about each plant used in the program. The plant objects are initialized with their
    scientific name and their code and then have different information added later.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.name: scientific name of the plant
    self.code: KEGG code for the plant
    self.genes: the gene entries for the plant
    self.ec_nums: the EC numbers parsed from the plants gene entries
    self.flavonoids: the list of flavonoids that the plant could potentially produce
    self.ec_counts: list of objects that hold the number of times each EC number appears
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    has_ec_count: used to determine whether or not a specific EC number is already in the list of EC counts
    incr_ec_count: used to increase the count for the EC count objects.
    """

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


class PathGene:
    """
    This object is used to hold Gene objects in a way such that they are sorted by the pathway from which they were
    found.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.path: the pathway that resulted in the gene entry
    self.genes: the list of gene entries from this pathway
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    """

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


class Gene:
    """
    This object holds data gathered from KEGG for each plant's pathway (like aip00491).
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.gene_id: the ID of the gene from a plant
    self.plant: the scientific name of the plant that has this gene
    self.plant_code: the KEGG code for the plant
    self.compound: the compound name listed in the entry
    self.ec_nums: the list of EC numbers found in the entry
    self.ortho: the KEGG orthology code for the compound
    self.path: the pathway where the gene was found
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    simple: returns a formatted string that contains information from the object
    no_plant: same as simple, but without including the plant name
    """

    def __init__(self, gene_id=None, plant=None, compound=None, ec_nums=None, ortho=None, path=None, plant_code=None):
        self.gene_id = gene_id if gene_id is not None else ' '
        self.plant = plant if plant is not None else ' '
        self.plant_code = plant_code if plant_code is not None else ' '
        self.compound = compound if compound is not None else ' '
        self.ec_nums = ec_nums if ec_nums is not None else []
        self.ortho = ortho if ortho is not None else ' '
        self.path = path if path is not None else ' '

    def __eq__(self, other):
        return self.gene_id == other.gene_id and \
               self.plant_code == other.plant_code and \
               self.ec_nums == other.ec_nums

    def is_in(self, items):
        for item in items:
            if self == item: return True
        return False

    def simple(self):
        return self.plant + ', ' + self.gene_id + ', ' + self.compound + ', ' + str(self.ec_nums) + ', ' + self.ortho

    def no_plant(self):
        return self.gene_id + ', ' + self.compound + ', ' + self.ec_nums + ', ' + self.ortho


class EcFastaCollection:
    """
    This object is used to hold the associated FASTA entries for any given EC number.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.ec_name: the EC number & name used when writing the file
    self.ec_entries: the list of associated FASTA entries (FastaEcEntry objects)
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    """

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
    """
    This object is a property of the Plant object and is used to hold each EC number and the number of times it occurs
    in gene entries of a given plant.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.number: the EC number
    self.count: number of times that the EC number shows up in gene entries.
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    """

    def __init__(self, number=None, count=None):
        self.number = number if number is not None else ' '
        self.count = count if number is not None else 0


class FastaEcEntry:
    """
    This object is a property of EcFastaCollection and contains the information for a specific FASTA entry.
    --------------------------------------------------------------------------------------------------------------------
    ATTRIBUTES
    self.gene_id: the gene ID associated with the sequence
    self.plant: the plant that the gene is from
    self.dna_seq: the dna sequence/FASTA entry for the specific gene
    --------------------------------------------------------------------------------------------------------------------
    FUNCTIONS
    __init__: constructor for the object
    __eq__: defines equality of the object
    is_in: determines if an identical or nearly identical object is already in the list
    simple: returns a formatted string
    """

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


data_lists = [ChemData(AGI, [], FN_AGI), ChemData(BUN, [], FN_BUN), ChemData(KXN, [], FN_KXN),
              ChemData(HWB, [], FN_HWB),
              ChemData(EC, [], FN_EC), ChemData(EGT, [], FN_EGT), ChemData(ERD, [], FN_ERD), ChemData(GC, [], FN_GC),
              ChemData(GEN, [], FN_GEN), ChemData(KMP, [], FN_KMP), ChemData(LU2, [], FN_LU2),
              ChemData(MYC, [], FN_MYC), ChemData(NAR, [], FN_NAR), ChemData(QUE, [], FN_QUER),
              ChemData(HCC, [], FN_HCC)]

import os
from parseutil import *

SEP = os.sep

FLAV_NAMES = 'pajson' + SEP + 'flav_names.json'
FLAV_REL = 'pajson' + SEP + 'flav_related.json'
FLAV_SYNS = 'pajson' + SEP + 'flav_syns.json'
LIG_ID_DICT = 'pajson' + SEP + 'lig_identifiers.json'
LIG_ID_PART = 'pajson' + SEP + 'lig_identifiers_part.json'
NCOMP_CODES = 'pajson' + SEP + 'chem_codes.json'
NCOMP_DICT = 'pajson' + SEP + 'chem_codes_names.json'
NPLANT_CODES = 'pajson' + SEP + 'organism_codes.json'
NPLANT_DICT = 'pajson' + SEP + 'organism_codes_names.json'
PDB_KNAP_IDS = 'pajson' + SEP + 'pdb_knap_ids.json'
PLANT_NAMES_CODES = 'pajson' + SEP + 'plant_names_codes.json'
PUB_PDB = 'pajson' + SEP + 'pub_pdb.json'
KS_FLAVS = 'pajson' + SEP + 'knapsack_flavs.json'

flav_list = get_json_data(FLAV_NAMES)
try:
    flav_relatives = get_json_data(FLAV_REL)
except UnicodeDecodeError:
    pass
flav_synonyms = get_json_data(FLAV_SYNS)
plant_dict = get_json_data(PLANT_NAMES_CODES)
flav_dict = get_json_data(KS_FLAVS)

CMP_CODES = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_codes.txt'
CMP_DICT = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_dict.txt'
CMP_IN = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_info.txt'
ORG_CODES = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'plant_codes.txt'
ORG_IN = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'npass_species.txt'
ORG_OUT = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'parsed_species.txt'
PAIRS_IN = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'plants_compounds.txt'
SELECTED = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'selected.txt'

RE_KS_PLANT = r'<td\sclass=org>\w+<.td><td\sclass=org>\w+<.td><td.class=org2>([\s\w]*)<'

RE_AROMA = r'<tr.id=.c[a-z]+Ar.+A.+C.+t.>\s*<th>A.+\sB.+\sC.+t</th>\s*<td>(\d+)</td>\s*</tr>'
RE_ATOM_COUNT = r'\s*<tr id=.ch.*A.*C.*t.>\s*<th>A.*C.*t</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=.ch.*Ch.*A.*C.*.>\s*<th>'
RE_BOND = r'<tr.id=.c[a-z]+B.+C.+t.>\s*<th>B.+\sC.+t</th>\s*<td>(\d*)</td>\s*</tr>\s*<tr.id=.ch.+A.+A.+C.+t.>'
RE_CHARGE = r'<tr id=.ch.+Fo.+C.+e.>\s*<th>Fo.+\sC.*e</th>\s*<td>([\w\.\-\+]+)</td>'
RE_CHEBI = r">CHEBI:([0-9]*)"
RE_CHEMBL = r">CHEMBL([0-9]*)"
RE_CHIRAL = r'<tr.id=.c[a-z]+C.+A.+t.>\s*<th>C[a-z]+\sAt.+Co.+t</th>\s*<td>(\d+)</td>\s*</tr>\s*<tr.id=.c.+B\w+C\w+t.>'
RE_COMP_ID = r'(NPC\S*)'
RE_COMP_NAME = r'NPC\S*\s*([^\t\r\n]*)'
RE_FORMULA = r'.*<tr id=.ch.*F.*a.>\s*<th>F.*a</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=.ch.*Mo.*W.*t\">'
RE_HAS_CAS = r'<tr><td class="d1"><a href=information\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td><td ' \
             r'class="d1">([0-9]*-[0-9]*-[0-9]*)<\/td><td class="d1">([^<]*)<\/td><td class="d1"> '
RE_INCHI = r"<tr id=\"chemicalInChI\"><th>InChI<\/th><td style=\"word-wrap: break-word\">(InChI=.*)</td></tr><tr " \
           r"id=\"chemicalInChIKey\"> "
RE_INCHI_KEY = r"<th>InChIKey<\/th><td>(.*)</td></tr></table></div><div class=\"col-md-4"
RE_KNAP_ENTRY = r'<tr>[.|\s]*<td\sc.*d1\">[.|\s]*<a.*blank\">(C[0-9]*)</a>[.|\s]*</td>[.|\s]*<td.*d1\">([ ' \
                r'0-9|-]*)</td>\s*<td\s.*d1\">(.*)</td>[.|\s]*<td\s.*d1\">([0-9|A-Z]+)</td>[.|\s]*<td\s.*d1\">([' \
                r'\d|\.]+).*(<[|/]td.*<[|/]tr>) '
RE_KNAP_NAME = r'<th class=\"inf\">Name</th>\s*<td colspan=\"4\" class=\"inf\">(.*)</td>\s*</tr>'
RE_KNAP_ORG = r'.*</td><td class=\"?org2\"?>([a-zA-Z0-9!@#$&()\\-`.+,/\"\s]*)'
RE_NAME_LONG = r'id=.ch[a-z]+I[a-z]+s.>\s*<th>Id.+rs</th>\s*<td.s.+=.word.+word.>(.*)</td></tr><tr.id=.ch[a-z]+Formula'
RE_NAME_REG = r'</ul>\s*</div>\s*<h1.id=.m.+eId.>.{3}</h1>\s*<h4 id=.m.+e.>\s*([^><!/].*)</h4>.*<div.class=.form-group'
RE_NO_CAS = r'<tr><td class="d1"><a h.*n\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td><td ' \
            r'class="d1"><\/td><td class="d1">([^<]*)<\/td><td class="d1"> '
RE_NUM_KNAP_RESULTS = r'Number of matched data :([0-9]*).*<br>'
RE_ORG_COMP_ID = r'-(NPC\S*)'
RE_ORG_ID = r'(NPO\S*)'
RE_ORG_NAME = r'NPO[^ \t\r]*\s*(\S*\s*\S*)\s*\bSpecies'
RE_ORG_NAME_PAIRS = r'(NPO\S*)\b-'
RE_PUBCHEM = r"https:\/\/pubchem\.ncbi\.nlm\.nih\.gov\/compound\/([0-9]*)"
RE_SMILES = r'id=.ch.*Is.*c.>\s*<th>Is.*S.*S</th>\s*<td.*\">(.*)</td>\s*</tr>\s*<tr.id=.ch.*I.*I.>\s*'
RE_WEIGHT = r'.*<tr id=.ch.*Mo.*W.*t.>\s*<th>Mo.*W.*t</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=\"ch.*lT.*e\">'

URL_KNAP_ORG = "http://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word="  # knapsack partial URL
URL_KNAP_CHEM = 'http://www.knapsackfamily.com/knapsack_core/information.php?sname=C_ID&word='
URL_PDB = "http://www.rcsb.org/ligand/"
URL_DBGET = 'https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+'

K_CHEB = 'CHEBI'
K_CHEM = 'CHEMBL'
K_FORM = 'FORMULA'
K_IKEY = 'INCHI_KEY'
K_INCHI = 'INCHI'
K_ISO = 'ISO_SMILES'
K_LONG = 'LONG'
K_NAME = 'COMMON'
K_PUB = 'PUBCHEM'
K_SMILE = 'SMILES'
S_IKEY = 'INCHIKEY'
S_INCHI = 'INCHI_CD'
S_METAB = 'metabolite'

ERR_VAL = 'ERR_WITH_REQUEST'
OUT_JSON = '/json'
OUT_TXT = '/txt'
PUBCHEM_PREFIX = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
SYN_HEAD = 'Synonyms (!! Separated)'
PROP_PAIRS = [('IUPACName', 'IUPAC'), ('InChI', 'InChI'), ('InChIKey', 'InChIKey'), ('CanonicalSMILES', 'CannonSMILES'),
              ('IsomericSMILES', 'IsoSMILES'), ('MolecularFormula', 'Formula'), ('MolecularWeight', 'Weight'),
              ('Charge', 'Charge'), ('Fingerprint2D', 'Fingerprint2D'), ('RotatableBondCount', '#RotatableBonds'),
              ('HeavyAtomCount', '#HeavyAtoms'), ('IsotopeAtomCount', '#IsotopeAtoms'),
              ('AtomStereoCount', '#StereoAtoms'), ('DefinedAtomStereoCount', '#DefStereoAtoms'),
              ('UndefinedAtomStereoCount', '#UndefStereoAtoms'), ('BondStereoCount', '#StereoBonds'),
              ('DefinedBondStereoCount', '#DefStereoBonds'), ('UndefinedBondStereoCount', '#UndefStereoBonds'),
              ('CovalentUnitCount', '#CovalentUnits'), ('HBondDonorCount', '#HBondDonors'),
              ('HBondAcceptorCount', '#HBondAcceptors')]
XREF_PAIRS = [('RN', 'CAS ID(S)'), ('MMDBID', 'MMDBID(S)'), ('ProteinGI', 'ProtGI(S)'), ('NucleotideGI', 'NTGI(S)'),
              ('TaxonomyID', 'TaxID(S)'), ('MIMID', 'MIMID(S)'), ('GeneID', 'GeneID(S)'), ('ProbeID', 'ProbeID(S)')]
MISC_PAIRS = [('ZINC', 'ZINC ID'), ('CTK', 'ChemTik ID'), ('CHEMBL', 'CHEMBL ID'), ('CHEBI:', 'CHEBI ID'),
              ('BDBM', 'BindingDB ID'), ('SCHEMBL', 'SureChEMBL ID'), ('DTXSID', 'EPA DSSTox')]

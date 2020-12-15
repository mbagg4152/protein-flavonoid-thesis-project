import urllib.request as url_req
import urllib.error as url_err
from util import *
from pathlib import Path
import os

ERR_VAL = 'ERR_WITH_REQUEST'
OUT_JSON = '/json'
OUT_TXT = '/txt'
PUBCHEM_PREFIX = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
SYN_HEAD = 'Synonyms ({} Separated)'
PROP_PAIRS = [('IUPACName', 'IUPAC'), ('InChI', 'InChI'), ('InChIKey', 'InChIKey'), ('CanonicalSMILES', 'CannonSMILES'),
              ('IsomericSMILES', 'IsoSMILES'), ('MolecularFormula', 'Formula'), ('MolecularWeight', 'Weight'), ('Charge', 'Charge'),
              ('Fingerprint2D', 'Fingerprint2D'), ('RotatableBondCount', '#RotatableBonds'), ('HeavyAtomCount', '#HeavyAtoms'),
              ('IsotopeAtomCount', '#IsotopeAtoms'), ('AtomStereoCount', '#StereoAtoms'), ('DefinedAtomStereoCount', '#DefStereoAtoms'),
              ('UndefinedAtomStereoCount', '#UndefStereoAtoms'), ('BondStereoCount', '#StereoBonds'),
              ('DefinedBondStereoCount', '#DefStereoBonds'), ('UndefinedBondStereoCount', '#UndefStereoBonds'),
              ('CovalentUnitCount', '#CovalentUnits'), ('HBondDonorCount', '#HBondDonors'), ('HBondAcceptorCount', '#HBondAcceptors')]
XREF_PAIRS = [('RN', 'CAS ID(S)'), ('MMDBID', 'MMDBID(S)'), ('ProteinGI', 'ProtGI(S)'), ('NucleotideGI', 'NTGI(S)'),
              ('TaxonomyID', 'TaxID(S)'), ('MIMID', 'MIMID(S)'), ('GeneID', 'GeneID(S)'), ('ProbeID', 'ProbeID(S)')]
MISC_PAIRS = [('ZINC', 'ZINC ID'), ('CTK', 'ChemTik ID'), ('CHEMBL', 'CHEMBL ID'), ('CHEBI:', 'CHEBI ID'), ('BDBM', 'BindingDB ID'),
              ('SCHEMBL', 'SureChEMBL ID'), ('DTXSID', 'EPA DSSTox')]

PROPS = [i[0] for i in PROP_PAIRS]
PROP_LBLS = [i[1] for i in PROP_PAIRS]
XREFS = [i[0] for i in XREF_PAIRS]
XREF_LBLS = [i[1] for i in XREF_PAIRS]
MISC = [i[0] for i in MISC_PAIRS]
MISC_LBLS = [i[1] for i in MISC_PAIRS]

ACT_PROP = '/property/' + ','.join(PROPS)
ACT_XREF = '/xrefs/' + ','.join(XREFS)
ACT_SYN = '/synonyms/'

pubchem_pdb = {}
total_out = 'PDB\tPC\t' + '\t'.join(PROP_LBLS) + '\t' + '\t'.join(XREF_LBLS) + '\t' + SYN_HEAD + '\t' + '\t'.join(MISC_LBLS) + '\n'
new_dir = '..' + SEP + 'misc_files' + SEP + 'pubchem' + SEP
images = new_dir + 'ligand_images' + SEP
prop_dir = new_dir + 'properties' + SEP
syn_dir = new_dir + 'synonyms' + SEP
xref_dir = new_dir + 'xref' + SEP
results = new_dir + 'results.tab'
Path(images).mkdir(parents=True, exist_ok=True)
Path(new_dir).mkdir(parents=True, exist_ok=True)
Path(prop_dir).mkdir(parents=True, exist_ok=True)
Path(syn_dir).mkdir(parents=True, exist_ok=True)
Path(xref_dir).mkdir(parents=True, exist_ok=True)

def main():
    global pubchem_pdb
    pubchem_pdb = get_json_data(FN_PUBCHEM_PDB)
    fetch_pubchem_data()
    with open(results, 'w') as file: file.write(total_out)

def fetch_pubchem_data():
    global total_out
    #####################################################################
    # Vals that need multiple columns
    # GeneID
    #####################################################################
    for key in pubchem_pdb:
        pc_id = pubchem_pdb.get(key)
        out = key + '\t' + pc_id + '\t'
        data = make_request(make_url(pc_id, ACT_PROP, OUT_JSON), 'PropertyTable', 'Properties', key, pc_id, 'PROPERTIES', prop_dir)
        if data != ERR_VAL:
            for prop in PROPS:
                val = data.get(prop)
                data_str = get_data_str(val, ',')
                out += data_str
                if len(data_str) > 5000: print('?????len of {} greater than 5K char for {}, {}'.format(prop, key, pc_id))
        else: out += ''.join(['NONE\t' for prop in PROPS])

        data = make_request(make_url(pc_id, ACT_XREF, OUT_JSON), 'InformationList', 'Information', key, pc_id, 'XREFS', xref_dir)
        if data != ERR_VAL:
            for xref in XREFS:
                val = data.get(xref)
                data_str = get_data_str(val, ',')
                out += data_str
                if len(data_str) > 5000: print('?????len of {} greater than 5K char for {}, {}'.format(xref, key, pc_id))
        else: out += ''.join(['NONE\t' for xref in XREFS])

        data = make_request(make_url(pc_id, ACT_SYN, OUT_JSON), 'InformationList', 'Information', key, pc_id, 'SYNONYMS', syn_dir)
        if data != ERR_VAL:
            val = data.get('Synonym')
            names = get_data_str(val, '{}')
            if len(names) > 5000: print('?????len of synonyms greater than 5K char for {}, {}'.format(key, pc_id))
            # out += names + '\t'
            out += ' ' + '\t'
            name_list = names.split('{}')
            for misc in MISC:
                tmp_list = []
                for name in name_list:
                    if name.strip()[0:len(misc)] == misc: tmp_list.append(name)
                if len(tmp_list): out += ','.join(tmp_list) + '\t'
                else: out += 'NONE\t'
        else: out += 'NONE\t' + ''.join(['NONE\t' for misc in MISC])

        total_out += out + '\n'

        if not os.path.exists(images + key + '.png'):
            try:
                url_req.urlretrieve(PUBCHEM_PREFIX + pc_id + '/PNG', images + key + '.png')
                print('got image for {} with pubchem id {}'.format(key, pc_id))
            except (url_err.HTTPError, url_err.URLError) as e:  request_err_handler('IMAGE', key, pc_id, e)

def make_request(url, key1, key2, pdb, pub, search, out_dir):
    has_file = False
    if os.path.exists(out_dir + pdb + '.json'):
        try:
            with open(out_dir + pdb + '.json') as tmp_json: content = json.load(tmp_json)
            has_file = True
        except FileNotFoundError as e: return request_err_handler(search, pdb, pub, e)
    else:
        try:
            content = url_req.urlopen(url).read().decode('utf-8').replace('\n', '')
            with open(out_dir + pdb + '.json', 'w') as tmp_file: tmp_file.write(content)
        except (url_err.HTTPError, url_err.URLError, ValueError, UnicodeDecodeError)as e:
            return request_err_handler(search, pdb, pub, e)

    try:
        if has_file:
            print('[data is from local file] got {} for {} with pubchem id {}'.format(search.lower(), pdb, pub))
            return content[key1][key2][0]
        else:
            print('got {} for {} with pubchem id {}'.format(search.lower(), pdb, pub))
            return json.loads(content)[key1][key2][0]
    except (IndexError, json.JSONDecodeError) as e: return request_err_handler(search, pdb, pub, e)

def request_err_handler(search, pdb, pub, err):
    print('!!!!NO {} FOR {} with PUBCHEM ID {}. Got {} with message {}'.format(search, pdb, pub, type(err).__name__, str(err)))
    return ERR_VAL

def get_data_str(val, sep):
    if isinstance(val, list):
        for i in range(0, len(val)):
            tmp = val[i]
            if not isinstance(tmp, str): val[i] = str(tmp)
        return sep.join(val) + '\t'
    else: return '{}'.format(val) + '\t'

def make_url(pc_id, action, output_form): return PUBCHEM_PREFIX + pc_id + action + output_form

if __name__ == '__main__': main()

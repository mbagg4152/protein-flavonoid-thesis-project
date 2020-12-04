from pathstrings import *
from miscvals import *
from util import *
import re
from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError

smiles_list = {}
re_chebi = r">CHEBI:([0-9]*)"
re_chembl = r">CHEMBL([0-9]*)"
re_inchi = r"<tr id=\"chemicalInChI\"><th>InChI<\/th><td style=\"word-wrap: break-word\">(InChI=.*)</td></tr>" \
           r"<tr id=\"chemicalInChIKey\">"
re_inchi_key = r"<th>InChIKey<\/th><td>(.*)</td></tr></table></div><div class=\"col-md-4"
re_knap_data = r"<tr class=\" row.*\"><td class=\"d1\"><a href=\"http://www\.knapsackfamily\.com/knapsack_core/" \
               r"information\.php\?word=(C[0-9]*)\" target=\"_blank\">C[0-9]*</a></td><td class=\"d1\">([\d|-]*)" \
               r"</td><td class=\"d1\">(.*)</td><td class=\"d1\">[a-zA-Z|\d]+</td><td class=\"d1\">[0-9|\.]*</td>" \
               r"<td class=\"d1\">"
re_pubchem_id = r"https:\/\/pubchem\.ncbi\.nlm\.nih\.gov\/compound\/([0-9]*)"
re_smiles = r"<tr id=\"chemicalIsomeric\"><th>Isomeric SMILES<\/th><td style=\"word-wrap: break-word\">(.*)</td>" \
            r"</tr><tr id=\"chemicalInChI\"><th>"

pdb_url = "http://www.rcsb.org/ligand/"
knap_url = "http://www.knapsackfamily.com/knapsack_core/result.php?sname=SMILES&word="

new_dir = '..' + SEP + 'misc_files' + SEP + 'smiles'
pages = new_dir + SEP + 'html_sites' + SEP
knap_data = new_dir + SEP + 'knapsack' + SEP
chem_info = new_dir + SEP + "chem_info.csv"
Path(new_dir).mkdir(parents=True, exist_ok=True)
Path(pages).mkdir(parents=True, exist_ok=True)
Path(knap_data).mkdir(parents=True, exist_ok=True)
def main():
    global smiles_list
    smiles_list = get_json_data(FN_SMILES)
    # get_extra_info()
    get_organism_data()

def get_organism_data():
    for key in smiles_list:
        print(smiles_list.get(key).get('SMILES'))
        item = smiles_list.get(key).get('SMILES').strip().replace(' ', '%20')
        tmp_url = knap_url + str(item)
        name = knap_data + key + ".txt"
        if not os.path.exists(name):
            try: request.urlretrieve(tmp_url, name)
            except HTTPError or URLError or http.client.InvalidURL or UnicodeEncodeError:
                print('err getting page')


def get_extra_info():
    out_str = ""
    for key in smiles_list:
        tmp_url = pdb_url + key
        name = pages + key + ".txt"

        if not os.path.exists(name):
            try: request.urlretrieve(tmp_url, name)
            except HTTPError or URLError as e:
                print('err getting page')
    for key in smiles_list:
        name = pages + key + ".txt"
        try: file = open(name, "r")
        except FileNotFoundError:
            print('no file')
            continue
        data = file.read()
        inchi, inchi_key, smiles, pubchem, chebi, chembl = 'NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE'
        try:
            inchi, inchi_key = re.findall(re_inchi, data)[0], re.findall(re_inchi_key, data)[0]
            smiles, pubchem = re.findall(re_smiles, data)[0], re.findall(re_pubchem_id, data)[0]
            chebi, chembl = re.findall(re_chebi, data)[0], re.findall(re_chembl, data)[0]
        except IndexError: pass
        out_str += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, smiles, inchi, inchi_key, pubchem, chebi, chembl)
        file.close()

    info_file = open(chem_info, 'w')
    info_file.write(out_str)
    info_file.close()

if __name__ == '__main__':
    main()

from http.client import InvalidURL
from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
from util import *
import re

smiles_list = {}
re_chebi = r">CHEBI:([0-9]*)"
re_chembl = r">CHEMBL([0-9]*)"
re_inchi = r"<tr id=\"chemicalInChI\"><th>InChI<\/th><td style=\"word-wrap: break-word\">(InChI=.*)</td></tr>" \
           r"<tr id=\"chemicalInChIKey\">"
re_inchi_key = r"<th>InChIKey<\/th><td>(.*)</td></tr></table></div><div class=\"col-md-4"
re_knap_data = r'<tr>[.|\s]*<td\sc.*d1\">[.|\s]*<a.*blank\">(C[0-9]*)</a>[.|\s]*</td>[.|\s]*<td.*d1\">([' \
               r'0-9|-]*)</td>\s*<td\s.*d1\">(.*)</td>[.|\s]*<td\s.*d1\">([0-9|A-Z]+)</td>[.|\s]*<td\s.*d1\">([' \
               r'\d|\.]+).*(<[|/]td.*<[|/]tr>)'
re_pubchem_id = r"https:\/\/pubchem\.ncbi\.nlm\.nih\.gov\/compound\/([0-9]*)"
re_smiles = r"<tr id=\"chemicalIsomeric\"><th>Isomeric SMILES<\/th><td style=\"word-wrap: break-word\">(.*)</td>" \
            r"</tr><tr id=\"chemicalInChI\"><th>"

re_knap_results = r'Number of matched data :([0-9]*).*<br>'

pdb_url = "http://www.rcsb.org/ligand/"

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
    get_organism_data()


def get_organism_data():
    for key in smiles_list:
        json_keys = ['SMILES', 'Name', 'Long', 'ISMILES', 'INCHI', 'INCHI_KEY']
        snames = ['SMILES', 'metabolite', 'metabolite', 'SMILES', 'INCHI_CD', 'INCHIKEY']
        words = []
        for sname in json_keys: words.append(smiles_list.get(key).get(sname).strip().replace(' ', '%20'))

        for i in range(0, len(json_keys)):
            encoded_words = u''.join(words[i]).encode('utf-8').strip()
            tmp_url = get_knap_url(snames[i], encoded_words)
            name = knap_data + key + ".html"
            if not os.path.exists(name):
                try: request.urlretrieve(tmp_url, name)
                except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
                    print('err getting page')
                    continue

            if not get_from_knapsack(tmp_url, name, key, json_keys[i]):
                if i != len(json_keys) - 1: os.remove(name)
                else:
                    file = open(name, 'w')
                    file.write('NO RESULTS')
                    file.close()
            else:
                break


def get_from_knapsack(tmp_url, name, key, prop):
    if not os.path.exists(name):
        try: request.urlretrieve(tmp_url, name)
        except HTTPError or URLError or InvalidURL or UnicodeEncodeError: print('err getting page')
    try:
        file = open(name, 'r')
        data = file.read()
        data = re.sub(r'<font color=#\S{6}>', '', data)
        data = re.sub(r'</font>', '', data)
        try:
            num_res = re.findall(re_knap_results, data)[0]
            if int(num_res) < 1:
                if prop == 'INCHI_KEY': print(
                    'no results for {} from value {}, saving data-less file.'.format(key, prop))
                else: print('no results for {} from value {}, deleting file.'.format(key, prop))
                return False
            else:
                out = ''
                short = False
                values = re.findall(re_knap_data, data)

                for item in values:
                    names = ' || '.join(item[2].split('<br>'))
                    last = item[5].replace('</td>', '').replace('</tr>', '')
                    out += '{} {} {} {} {} {}\n'.format(item[0], item[1], names, item[3], item[4], str(last))
                print('for {}\n{}'.format(key, out))
        except IndexError: return False
        file.close()
    except FileNotFoundError:
        print('no file for {} from value {}'.format(key, prop))
        return False

    return True


def get_knap_url(sname, word):
    return "http://www.knapsackfamily.com/knapsack_core/result.php?sname={}&word={}".format(sname, word)


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

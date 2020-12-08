from http.client import InvalidURL
from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
from util import *
import re
import glob

json_keys = ['SMILES', 'Name', 'Long', 'ISMILES', 'INCHI', 'INCHI_KEY']
snames = ['SMILES', 'metabolite', 'metabolite', 'SMILES', 'INCHI_CD', 'INCHIKEY']
smiles_list = {}
knap_ids = []
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
re_knap_name = r'<th class=\"inf\">Name</th>\s*<td colspan=\"4\" class=\"inf\">(.*)</td>\s*</tr>'
re_knap_results = r'Number of matched data :([0-9]*).*<br>'
re_knap_org = r'.*</td><td class=\"?org2\"?>([a-zA-Z0-9!@#$&()\\-`.+,/\"\s]*)'
pdb_url = "http://www.rcsb.org/ligand/"
re_search_type = r'.*input?\stype?\s=?\s<font?\sclass=\"iw\">\s*(\S*)?\s,?\s</font>'
new_dir = '..' + SEP + 'misc_files' + SEP + 'smiles'
pages = new_dir + SEP + 'html_sites' + SEP
match_dir = new_dir + SEP + 'matches' + SEP
knap_dir = new_dir + SEP + 'knap' + SEP
knap_data = {}
matches = ''
chem_info = new_dir + SEP + "chem_info.csv"
match_info = new_dir + SEP + 'matches.csv'
knap_id_info = new_dir + SEP + 'knap_ids.csv'
knap_orgs = new_dir + SEP + 'knap_orgs.csv'
Path(new_dir).mkdir(parents=True, exist_ok=True)
Path(pages).mkdir(parents=True, exist_ok=True)
Path(match_dir).mkdir(parents=True, exist_ok=True)
Path(knap_dir).mkdir(parents=True, exist_ok=True)

def main():
    global smiles_list, knap_ids
    smiles_list = get_json_data(FN_SMILES)
    # get_organism_data()
    knap_ids = get_json_data(FN_KNAP)

    with open(knap_id_info, 'w') as outfile:
        out = ''
        for knap in knap_ids: out += knap + '\n'
        outfile.write(out)
    get_knap_data()

def is_ascii(s):
    return all(ord(c) < 128 for c in s)


def get_organism_data():
    for key in smiles_list:
        words = []
        for sname in json_keys: words.append(smiles_list.get(key).get(sname).strip().replace(' ', '%20'))

        for i in range(0, len(json_keys)):
            if not is_ascii(words[i]): words[i] = ''.join(str(ord(c)) for c in words[i])

            encoded_words = u''.join(words[i]).encode('utf-8').strip().decode('utf-8')
            tmp_url = get_knap_url(snames[i], encoded_words)
            name = match_dir + key + ".html"
            files = ''.join(glob.glob(match_dir + key + '_*html')).strip()
            if '_NONE' in files:
                break

            if not (os.path.exists(name)):
                print('NEED TO GET FILE {}'.format(files))
                if not (os.path.exists(name)):
                    try:
                        request.urlretrieve(tmp_url, name)
                        print('{} {}'.format(tmp_url, name))
                    except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
                        print('err getting page')
                        return False
            else:
                if not has_codes_from_knapsack(name, key, json_keys[i], tmp_url):
                    if snames[i] != 'INCHIKEY' and json_keys[i] != 'INCHI_KEY':
                        try: os.remove(name)
                        except FileNotFoundError:
                            print('##NO DELETE')
                    else:
                        new_name = match_dir + key + '_NONE' + ".html"
                        has_codes_from_knapsack(new_name, key, json_keys[i], tmp_url)
                        try: os.remove(name)
                        except FileNotFoundError:
                            print('##NO DELETE')

    with open(match_info, 'w') as match_file:
        match_file.write(matches)

def get_knap_data():
    global knap_ids
    out = ''
    for kid in knap_ids:
        tmp_url = 'http://www.knapsackfamily.com/knapsack_core/information.php?sname=C_ID&word=' + kid.strip()
        tmp_name = knap_dir + kid.strip() + '.html'
        if not (os.path.exists(tmp_name)):
            try: request.urlretrieve(tmp_url, tmp_name)
            except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
                print('err getting page')
                return False
        try:
            with open(tmp_name, 'r') as tmp_file:
                data = tmp_file.read()
                orgs = '\t'.join(re.findall(re_knap_org, data))
                name = kid
                if len(re.findall(re_knap_name, data)):
                    name = ' || '.join(re.findall(re_knap_name, data)[0].split('<br>')) + ':'

                out += name + '\t' + orgs + '\n'


        except FileNotFoundError:
            continue
    with open(knap_orgs, 'w') as file:
        file.write(out)

def has_codes_from_knapsack(name, key, prop, tmp_url):
    global matches, knap_data, knap_ids
    if not (os.path.exists(name)):
        try:
            request.urlretrieve(tmp_url, name)
            print('{} {}'.format(tmp_url, name))
        except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
            print('err getting page')
            return False
    try:
        with open(name, 'r') as file:
            data = file.read()
            data = re.sub(r'<font color=#\S{6}>', '', re.sub(r'</font>', '', data))
            try:
                num_res = re.findall(re_knap_results, data)[0]
                if int(num_res) < 1: return False
                else:
                    out = ''
                    values = re.findall(re_knap_data, data)

                    for item in values:
                        names = ' || '.join(item[2].split('<br>'))
                        last = item[5].replace('</td>', '').replace('</tr>', '')
                        matches += '{} {} {} {} {} {}\n'.format(item[0], item[1], names, item[3], item[4],
                                                                str(last))
                        if item[0] not in knap_ids: knap_ids.append(item[0])
                    return True
            except IndexError: return False


    except FileNotFoundError:
        print('no file for {} from property {}, {}'.format(key, prop, name))
        with open(name, 'w')as failed:
            failed.write('NOTHING')
        return False


def get_knap_url(sname, word):
    return "http://www.knapsackfamily.com/knapsack_core/result.php?sname={}&word=".format(sname) + word


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

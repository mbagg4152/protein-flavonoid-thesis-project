from http.client import InvalidURL
from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
# from util import *
from re import findall, sub
from parsersconstants import *
import glob

knap_ids = []
org_list = []
json_keys = [K_INCHI, K_IKEY, K_SMILE, K_ISO, K_LONG, K_NAME]
snames = [S_INCHI, S_IKEY, K_SMILE, K_SMILE, S_METAB, S_METAB]
knap_data = {}
part_list = {}
smiles_list = {}
matches = 'LigID\tC_ID\tMatch\tName\tLong\tKSNames\tForm\tKSForm\tSMILES\tIsomeric\tKS_SMILES\tINCHI\tINCHIKEY\t' \
          'CAS_ID\tCHEMBL\tCHEBI\tPUBCHEM'

new_dir = '..' + SEP + 'misc_files' + SEP + 'smiles'

knap_dir = new_dir + SEP + 'knap' + SEP
match_dir = new_dir + SEP + 'matches' + SEP
pages = new_dir + SEP + 'html_sites' + SEP
chem_info = new_dir + SEP + "chem_info.tab"
knap_id_info = new_dir + SEP + 'knap_ids.tab'
knap_orgs = new_dir + SEP + 'knap_orgs.tab'
match_info = new_dir + SEP + 'matches.tab'
org_names = new_dir + SEP + 'orgs.tab'
Path(knap_dir).mkdir(parents=True, exist_ok=True)
Path(match_dir).mkdir(parents=True, exist_ok=True)
Path(new_dir).mkdir(parents=True, exist_ok=True)
Path(pages).mkdir(parents=True, exist_ok=True)


def main():
    global smiles_list, knap_ids, part_list
    smiles_list = get_json_data(lig_id_dict)
    part_list = get_json_data(lig_id_part)
    knap_ids = get_json_data(pdb_knap_ids)
    get_knapsack_codes()
    with open(knap_id_info, 'w') as outfile:
        out = ''
        for knap in knap_ids: out += knap + '\n'
        outfile.write(out)
    # get_extra_info()
    with open(match_info, 'w') as outfile:
        outfile.write(matches)
    get_parse_knap_pages()


def get_knapsack_codes():
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
                print('no matches at all. {} has a _NONE file'.format(key))
                break

            if not (os.path.exists(name)):
                if not (os.path.exists(name)):
                    try:
                        request.urlretrieve(tmp_url, name)
                    except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
                        print('err getting page')
                        continue

            if not knap_code_helper(name, key, json_keys[i], tmp_url):
                if i != len(json_keys) - 1:
                    print('no match {} {}'.format(key, json_keys[i]))
                    try:
                        os.remove(name)
                    except FileNotFoundError:
                        print('!!!NO DELETE')
                else:
                    print('no matches at all {} {} ... making _NONE file.'.format(key, json_keys[i]))
                    new_name = match_dir + key + '_NONE' + ".html"
                    knap_code_helper(new_name, key, json_keys[i], tmp_url)
                    try:
                        os.remove(name)
                    except FileNotFoundError:
                        print('!!!NO DELETE')
            else:
                print('MATCH {} {}'.format(key, json_keys[i]))
                break

    with open(match_info, 'w') as match_file:
        match_file.write(matches)


def get_parse_knap_pages():
    global knap_ids, org_list
    out = ''
    for key in smiles_list:
        tmp_dict = smiles_list.get(key)
        kid = tmp_dict.get('KNAP_ID') or 'NONE'
        if kid == 'NONE': continue
        tmp_url = 'http://www.knapsackfamily.com/knapsack_core/information.php?sname=C_ID&word=' + kid.strip()
        tmp_name = knap_dir + kid.strip() + '.html'
        if not (os.path.exists(tmp_name)):
            try:
                request.urlretrieve(tmp_url, tmp_name)
            except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
                print('err getting page')
                return False
        try:
            with open(tmp_name, 'r') as tmp_file:
                data = tmp_file.read()
                found = findall(re_knap_org, data)
                orgs = '\t'.join(found)
                for org in found:
                    if org not in org_list: org_list.append(org)
                name = key
                if len(findall(re_knap_name, data)): name = ' || '.join(findall(re_knap_name, data)[0].split('<br>'))
                out += key + '\t' + kid + '\t' + name + '\t' + orgs + '\n'
        except FileNotFoundError:
            continue
    with open(knap_orgs, 'w') as file:
        file.write(out)
    with open(org_names, 'w') as file:
        org_out = ''
        for org in org_list: org_out += org + '\n'
        file.write(org_out)


def knap_code_helper(out_path, key, prop, tmp_url):
    global matches, knap_data, knap_ids
    if not (os.path.exists(out_path)):
        try:
            request.urlretrieve(tmp_url, out_path)
        except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
            print('err getting page')
            return False
    try:
        with open(out_path, 'r') as file:
            data = file.read()
            data = sub(r'<font color=#\S{6}>', '', sub(r'</font>', '', data))
            try:
                num_res = findall(re_num_knap_results, data)[0]
                if int(num_res) < 1:
                    return False
                else:
                    values = findall(re_knap_entry, data)
                    tmp_entry = smiles_list.get(key) or 'NONE'
                    name = tmp_entry.get(K_NAME) or 'NONE'
                    lname = tmp_entry.get(K_LONG) or 'NONE'
                    smiles = tmp_entry.get(K_SMILE) or 'NONE'
                    iso = tmp_entry.get(K_ISO) or 'NONE'
                    form = tmp_entry.get(K_FORM) or 'NONE'
                    inchi = tmp_entry.get(K_INCHI) or 'NONE'
                    ikey = tmp_entry.get(K_IKEY) or 'NONE'
                    pub = tmp_entry.get(K_PUB) or 'NONE'
                    chebi = tmp_entry.get(K_CHEB) or 'NONE'
                    chembl = tmp_entry.get(K_CHEM) or 'NONE'
                    weight = tmp_entry.get('WEIGHT') or 'NONE'
                    for item in values:
                        ks_smiles = sub(r'</td>', '', sub(r'</tr>', '', item[5])).strip()
                        if ks_smiles == '':
                            if prop == K_SMILE:
                                ks_smiles = ''.join(findall(r'.*word=(.*)', tmp_url))
                            else:
                                ks_smiles = 'NO SMILES LISTED DIRECTLY'
                        c_id = item[0]
                        cas_id = item[1]
                        ks_names = ' || '.join(item[2].split('<br>'))
                        ks_form = item[3]
                        ks_weight = item[4]
                        has_exact = False
                        if not has_exact:
                            split_names = ks_names.split(' || ')
                            for compound in split_names:
                                if bskin(compound) == bskin(name) or bskin(compound) == bskin(lname):
                                    has_exact = True
                                    break
                            if not has_exact:
                                has_exact = smiles.strip() == ks_smiles.strip() or iso.strip() == ks_smiles.strip()

                        if bskin(ks_form) != bskin(form):
                            has_exact = False
                        if has_exact:
                            # LigID C_ID Match Name Long KSNames Form KSForm SMILES Isomeric KS_SMILES INCHI INCHIKEY
                            # CAS_ID CHEMBL CHEBI PUBCHEM
                            # LigID, C_ID & Match added later
                            tmp_out = '\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}' \
                                      '\t{}\t{}\t{}\t{}\t{}'.format(name, lname, ks_names, form,
                                                                    ks_form, smiles, iso, ks_smiles, inchi, ikey,
                                                                    cas_id, chembl, chebi, pub)

                            if item[0] not in knap_ids:
                                matches += '\n' + key + '\t' + c_id + '\t' + prop + '\t' + tmp_out
                                knap_ids.append(item[0])

                    return True
            except IndexError:
                return False
    except FileNotFoundError:
        print('no file for {} from property {}, {}'.format(key, prop, out_path))
        with open(out_path, 'w')as failed:
            failed.write('NOTHING')
        return False


def get_knap_url(sname, word):
    return "http://www.knapsackfamily.com/knapsack_core/result.php?sname={}&word=".format(sname) + word


def get_extra_info():
    out_str = ""
    for key in smiles_list:
        tmp_url = PDB_URL + key
        name = pages + key + ".txt"
        if not os.path.exists(name):
            try:
                request.urlretrieve(tmp_url, name)
            except HTTPError or URLError as e:
                print('err getting page')
    for key in smiles_list:
        name = pages + key + ".txt"
        try:
            file = open(name, "r")
        except FileNotFoundError:
            print('no file')
            continue
        data = file.read()
        aroma = ''.join(findall(re_aroma, data) or ['NONE'])
        bond = ''.join(findall(re_bond, data) or ['NONE'])
        chebi = ''.join(findall(re_chebi, data)[0]) if len(findall(re_chebi, data)) else 'NONE'
        chembl = ''.join(findall(re_chembl, data)[0]) if len(findall(re_chembl, data)) else 'NONE'
        chiral = ''.join(findall(re_chiral, data) or ['NONE'])
        count = ''.join(findall(re_atom_count, data) or ['NONE'])
        f_charge = ''.join(findall(re_charge, data) or ['NONE'])
        form = sub(' ', '', ''.join(sub('<sub>', '', sub('</sub>', '', (''.join(findall(re_form, data) or ['NONE']))))))
        inchi = ''.join(findall(re_inchi, data) or ['NONE'])
        inchi_key = ''.join(findall(re_inchi_key, data) or ['NONE'])
        pubchem = ''.join(findall(re_pubchem, data)[0]) if len(findall(re_pubchem, data)) else 'NONE'
        smiles = ''.join(findall(re_smile, data)[0]) if len(findall(re_smile, data)) else 'NONE'
        weight = findall(re_weight, data)[0] if len(findall(re_weight, data)) else 'NONE'

        common = ''.join(findall(re_name, data)[0]) if len(findall(re_name, data)) else 'NONE'
        long = ''.join(findall(re_long, data)[0]) if len(findall(re_long, data)) else 'NONE'

        out_str += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}" \
                   "\t{}\t{}\t{}\n".format(key, smiles, inchi, inchi_key, pubchem, chebi, chembl,
                                           form, weight, f_charge, count, chiral, bond, aroma)

        file.close()

    info_file = open(chem_info, 'w')
    info_file.write(out_str)
    info_file.close()


def is_ascii(s):
    return all(ord(c) < 128 for c in s)


def bskin(line): return ''.join(line.strip().upper())


if __name__ == '__main__':
    main()

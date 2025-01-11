from urllib import request
from urllib.error import HTTPError, URLError

from paconstants import *

new_dir_plants = '..' + SEP + 'misc_files' + SEP + 'knapsack-plants'
comp_dir = '..' + SEP + 'misc_files' + SEP + 'knapsack-comps'
ks_plants_file = new_dir_plants + SEP + 'processed_plants.tsv'
ks_plant_all_comps_file = new_dir_plants + SEP + 'plants_all_comps.tsv'
ks_chem_file = comp_dir + SEP + 'processed_compounds.tsv'
plant_flavs = {}
all_plant_comps = {}
flav_orgs = {}
PLANT_PARSE = True

def main():
    try: os.mkdir(new_dir_plants)
    except FileExistsError: pass
    try: os.mkdir(new_dir_plants + SEP + 'raw-html')
    except FileExistsError: pass

    try: os.mkdir(comp_dir)
    except FileExistsError: pass
    try: os.mkdir(comp_dir + SEP + 'raw-html')
    except FileExistsError: pass

    # do_chem_parsing()
    do_plant_parsing()

def get_knap_file(name, url):
    if not os.path.exists(name):
        try:
            print('getting ' + url)
            request.urlretrieve(url, name)
        except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
            print('err getting page')

def do_chem_parsing():
    with open(ks_chem_file, 'w') as f:
        for key in knap_flav_dict:
            tmp_name = knap_flav_dict[key]
            tmp_url = URL_KNAP_CHEM + key.strip()
            tmp_fname = comp_dir + SEP + 'raw-html' + SEP + tmp_name + '.txt'
            get_knap_file(tmp_fname, tmp_url)
            parse_chem_file(tmp_fname, tmp_name)

        ks_str = ''

        for key in flav_orgs:
            if len(flav_orgs[key]) > 0:
                for item in flav_orgs[key]:
                    ks_str += item + '\t' + key + '\n'

        f.write(ks_str)  # write the formatted string to the file

def parse_chem_file(file_name, chem_name):
    file = open(file_name, 'r')
    lines = file.readlines()  #
    flav_orgs[chem_name] = []
    for line in lines:
        tmp_plant = re.findall(RE_KS_PLANT, line)
        if len(tmp_plant) > 0:
            flav_orgs[chem_name].append(tmp_plant[0])

def do_plant_parsing():
    global plant_flavs

    for key in plant_dict:
        tmp_val = plant_dict[key].split()[0] + ' ' + plant_dict[key].split()[1]

        # plant_flavs[tmp_val] = []
        tmp_url = URL_KNAP_ORG + tmp_val  # fill out URL for the plant
        tmp_url = tmp_url.replace(' ', '%20')
        tmp_file_name = new_dir_plants + SEP + 'raw-html' + SEP + key + ".txt"
        get_knap_file(tmp_file_name, tmp_url)  # download the page for the plant
        parse_plant_file(tmp_file_name, tmp_val)  # parse the downloaded file

    with open(ks_plants_file, 'w') as f:
        ks_str = ''
        for key in plant_flavs:
            if len(plant_flavs[key]) > 0:  # at least one entry was found for the plant
                ks_str += key + '\t' + '\t'.join(plant_flavs.get(key)) + '\n'
            else: ks_str += key + '\n'

        f.write(ks_str)  # write the formatted string to the file
    with open(ks_plant_all_comps_file, 'w') as f:
        ks_str = ''
        for key in all_plant_comps:
            if len(all_plant_comps[key]) > 0:  # at least one entry was found for the plant
                all_plant_comps[key].sort()
                for item in all_plant_comps.get(key):
                    ks_str += key + '\t' + item + '\n'
                # ks_str += key + '\t' + '\t'.join(all_plant_comps.get(key)) + '\n'
            else: ks_str += key + '\n'

        f.write(ks_str)  # write the formatted string to the file

def parse_plant_file(file_name, plant_name):
    global plant_flavs, all_plant_comps
    with open(file_name, 'r') as file:
        lines = file.readlines()  # each metabolite line begins with this HTML tag

        for line in lines:
            has_cas = re.findall(RE_HAS_CAS, line)
            no_cas = re.findall(RE_NO_CAS, line)

            if len(has_cas) > 0:
                tmp_id = has_cas[0][0]
                tmp_name = has_cas[0][2]
                knap_org_name = (has_cas[0][3] + '' + has_cas[0][4] + ' ' + has_cas[0][5]).strip()
            elif len(no_cas) > 0:
                tmp_id = no_cas[0][0]
                tmp_name = no_cas[0][1]
                knap_org_name = (no_cas[0][2] + '' + no_cas[0][3] + ' ' + no_cas[0][4]).strip()
            else: continue

            if knap_org_name not in plant_flavs.keys(): plant_flavs[knap_org_name] = []
            if knap_org_name not in all_plant_comps.keys(): all_plant_comps[knap_org_name] = []
            for key in knap_flav_dict:
                if tmp_name not in all_plant_comps.get(knap_org_name):
                    all_plant_comps[knap_org_name].append(tmp_name)
                if key.strip() in tmp_id.strip():
                    plant_flavs[knap_org_name].append(knap_flav_dict.get(key))

if __name__ == '__main__':
    main()

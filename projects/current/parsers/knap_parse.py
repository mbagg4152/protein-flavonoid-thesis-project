import sys
import os
import re
from paconstants import *
from urllib import request
from urllib.error import HTTPError, URLError
import urllib.parse

new_dir_plants = '..' + SEP + 'misc_files' + SEP + 'ks_plant_parse'
comp_dir = '..' + SEP + 'misc_files' + SEP + 'ks_comp_parse'
ks_plants_csv = new_dir_plants + SEP + '000_knapsack_data.csv'
ks_chem_csv = comp_dir + SEP + '000_ks_chem_data.csv'
plant_flavs = {}
flav_orgs = {}
PLANT_PARSE = True


def main():
    try: os.mkdir(new_dir_plants)
    except FileExistsError: pass
    try: os.mkdir(comp_dir)
    except FileExistsError: pass
    for key in flav_dict:
        tmp_name = flav_dict[key]
        tmp_url = URL_KNAP_CHEM + key.strip()
        tmp_fname = comp_dir + SEP + tmp_name + '.txt'
        get_knap_file(tmp_fname, tmp_url)
        parse_chem_file(tmp_fname, tmp_name)

    ks_file = open(ks_chem_csv, 'w')  # make file for outputting the knapsack data
    ks_str = ''

    for key in flav_orgs:
        if len(flav_orgs[key]) > 0:
            for item in flav_orgs[key]:
                ks_str += item + ',' + key + '\n'

    for key in plant_dict:
        tmp_val = plant_dict[key]
        tmp_url = URL_KNAP_ORG + tmp_val  # fill out URL for the plant
        tmp_url = tmp_url.replace(' ', '%20')
        tmp_file_name = new_dir_plants + SEP + key + ".txt"  # temporary filename for the downloaded page
        get_knap_file(tmp_file_name, tmp_url)  # download the page for the plant
        parse_plant_file(tmp_file_name, tmp_val)  # parse the downloaded file

    ks_file.write(ks_str)  # write the formatted string to the file
    ks_file.close()

    # ks_file = open(ks_plants_csv, 'w')  # make file for outputting the knapsack data
    # ks_str = ''
    #
    # for key in plant_flavs:
    #     if len(plant_flavs[key]) > 0:  # at least one entry was found for the plant
    #         for item in plant_flavs[key]:
    #             # add to the output string such that the data is written in CSV format
    #             ks_str += key + ',' + item[0] + ',\"' + item[1] + '\"\n'
    #     else:
    #         ks_str += key + '\n'
    #
    # ks_file.write(ks_str)  # write the formatted string to the file
    # ks_file.close()


def get_knap_file(name, url):
    if not os.path.exists(name):
        try:
            print('getting ' + url)
            request.urlretrieve(url, name)
        except HTTPError or URLError or InvalidURL or UnicodeEncodeError:
            print('err getting page')


def parse_chem_file(file_name, chem_name):
    file = open(file_name, 'r')
    lines = file.readlines()  #
    flav_orgs[chem_name] = []
    for line in lines:
        tmp_plant = re.findall(RE_KS_PLANT, line)
        if len(tmp_plant) > 0:
            flav_orgs[chem_name].append(tmp_plant[0])


def parse_plant_file(file_name, plant_name):
    file = open(file_name, 'r')
    lines = file.readlines()  # each metabolite line begins with this HTML tag
    plant_flavs[plant_name] = []
    for line in lines:
        tmp_line = line
        flav_match = 0
        syn_match = 0
        rel_match = 0
        # some entries don't have a CAS ID. if nothing is found using regular expression testing, then the
        # flavonoid name will follow right after the knapsack ID
        tmp_has_cas = re.findall(RE_HAS_CAS, line)
        tmp_no_cas = re.findall(RE_NO_CAS, line)
        tmp_list = []
        if len(tmp_no_cas) > 0:
            tmp_id = tmp_no_cas[0][0]
            tmp_name = tmp_no_cas[0][1]
        elif len(tmp_has_cas) > 0:
            tmp_id = tmp_has_cas[0][0]
            tmp_name = tmp_has_cas[0][2]
        else:
            continue
        tmp_list.append(tmp_id)
        tmp_list.append(tmp_name)
        for key in flav_list:
            tmp_cleaned = tmp_name.replace('(', '').replace(')', '').replace('+', '').replace('-', '')
            tmp_list[1] = tmp_cleaned
            if key.upper().strip() in tmp_cleaned.upper().strip():
                flav_match += 1
                plant_flavs[plant_name].append(tmp_list)
                continue

        if flav_match == 0:
            for key in flav_synonyms:
                for flav in flav_synonyms[key]:
                    if flav.upper().strip() == tmp_name.upper().strip():
                        syn_match += 1
                        tmp_list[1] = '{' + key + ' ALT NAME} ' + tmp_name
                        plant_flavs[plant_name].append(tmp_list)
                    elif flav in tmp_name:
                        tmp_list[1] = '{' + key + ' REL} ' + tmp_name
                        plant_flavs[plant_name].append(tmp_list)

            # if syn_match == 0:
            #     for key in flav_relatives:
            #         for flav in flav_relatives[key]:
            #             if flav.upper().strip() in tmp_name.upper().strip():
            #                 rel_match += 1
            #                 tmp_list[1] = '{' + key + ' REL} ' + tmp_name
            #                 plant_flavs[plant_name].append(tmp_list)
            #     if rel_match == 0:
            #         for flav in flav_list:
            #             if flav.upper().strip() in tmp_name.upper().strip():
            #                 tmp_list[1] = '{' + flav + ' IN NAME} ' + tmp_name
            #                 plant_flavs[plant_name].append(tmp_list)

    file.close()


if __name__ == '__main__':
    main()

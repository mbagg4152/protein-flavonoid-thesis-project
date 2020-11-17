from jsondata import *
from pathstrings import SEP
from miscvals import *
import os
import re

new_dir = '..' + SEP + 'misc_files' + SEP + 'knapsack_dir'
ks_data_name = new_dir + SEP + '000_knapsack_data.csv'
plant_flavs = {}


def main():
    try: os.mkdir(new_dir)
    except FileExistsError: pass

    # parsing file of plant codes and names
    for key in plant_dict_reg:
        tmp_val = plant_dict_reg[key]  # plant name
        tmp_url = URL_KNAP + tmp_val + "'"  # fill out URL for the plant
        tmp_file_name = new_dir + SEP + key + ".txt"  # temporary filename for the downloaded page
        wget_download(tmp_file_name, tmp_url)  # download the page for the plant
        parse_file(tmp_file_name, tmp_val)  # parse the downloaded file

    ks_file = open(ks_data_name, 'w')  # make file for outputting the knapsack data
    ks_str = ''
    for key in plant_flavs:
        if len(plant_flavs[key]) > 0:  # at least one entry was found for the plant
            for item in plant_flavs[key]:
                # add to the output string such that the data is written in CSV format
                ks_str += key + ',' + item[0] + ',\"' + item[1] + '\"\n'
        else: ks_str += key + '\n'
    ks_file.write(ks_str)  # write the formatted string to the file
    ks_file.close()


def wget_download(name, url):
    tmp_cmd = "wget -O " + name + ' ' + url  # custom wget command
    if not os.path.exists(name): os.system(tmp_cmd)  # use the system's wget command


def parse_file(file_name, plant_name):
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
        tmp_has_cas = re.findall(re_line_with_cas, line)
        tmp_no_cas = re.findall(re_line_no_cas, line)
        tmp_list = []
        if len(tmp_no_cas) > 0:
            tmp_id = tmp_no_cas[0][0]
            tmp_name = tmp_no_cas[0][1]
        elif len(tmp_has_cas) > 0:
            tmp_id = tmp_has_cas[0][0]
            tmp_name = tmp_has_cas[0][2]
        else: continue
        tmp_list.append(tmp_id)
        tmp_list.append(tmp_name)
        for key in flav_list:
            tmp_cleaned = tmp_name.replace('(', '').replace(')', '').replace('+', '').replace('-', '')
            tmp_list[1] = tmp_cleaned
            if key.upper().strip() == tmp_cleaned.upper().strip():
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

            if syn_match == 0:
                for key in flav_relatives:
                    for flav in flav_relatives[key]:
                        if flav.upper().strip() in tmp_name.upper().strip():
                            rel_match += 1
                            tmp_list[1] = '{' + key + ' REL} ' + tmp_name
                            plant_flavs[plant_name].append(tmp_list)
                if rel_match == 0:
                    for flav in flav_list:
                        if flav.upper().strip() in tmp_name.upper().strip():
                            tmp_list[1] = '{' + flav + ' IN NAME} ' + tmp_name
                            plant_flavs[plant_name].append(tmp_list)

    file.close()


if __name__ == '__main__':
    main()

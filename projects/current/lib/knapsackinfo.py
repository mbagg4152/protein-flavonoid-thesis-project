from jsondata import *
from pathstrings import SEP
from miscvals import *
import os
import re

new_dir = '..' + SEP + 'misc_files' + SEP + 'knapsack_dir'
ks_data_name = new_dir + SEP + 'knapsack_data.csv'
plant_flavs = {}
cas_id_form = r"([0-9]*-[0-9]*-[0-9]*)"  # regular expression for CAS ID
re_line_with_cas = r'<tr><td class="d1"><a href=information\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td>' \
                   r'<td class="d1">([0-9]*-[0-9]*-[0-9]*)<\/td><td class="d1">([^<]*)<\/td><td class="d1">'
re_line_no_cas = r'<tr><td class="d1"><a href=information\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td><td ' \
                 r'class="d1"><\/td><td class="d1">([^<]*)<\/td><td class="d1">'

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

    try:  # make output file if it does not exist
        tmp_file = open(ks_data_name, 'x')
        tmp_file.close()
    except FileExistsError or PermissionError: print('')

    ks_file = open(ks_data_name, 'w')  # make file for outputting the knapsack data
    ks_str = ''
    for key in plant_flavs:
        if len(plant_flavs[key]) > 0:  # at least one entry was found for the plant

            for item in plant_flavs[key]:
                # some entries don't have a CAS ID. if nothing is found using regular expression testing, then the
                # flavonoid name will follow right after the knapsack ID
                check = re.findall(cas_id_form, item[2])  # see if the line has a CAS ID

                # add to the output string such that the data is written in CSV format
                if len(check) < 1: ks_str += key + ',' + item[0] + ',\"' + item[2] + '\"\n'
                else: ks_str += key + ',' + item[0] + ',\"' + item[3] + '\"\n'
        else: ks_str += key + ''
        ks_str += '\n'
    ks_file.write(ks_str)  # write the formatted string to the file
    ks_file.close()


def wget_download(name, url):
    # make sure temp file exists
    try:
        tmp_file = open(name, 'x')
        tmp_file.close()
    except FileExistsError: pass
    tmp_cmd = "wget -O " + name + ' ' + url  # custom wget command
    os.system(tmp_cmd)  # use the system's wget command


def parse_file(file_name, plant_name):
    file = open(file_name, 'r')
    lines = file.read().split('<tr')  # each metabolite line begins with this HTML tag
    plant_flavs[plant_name] = []
    for line in lines:
        tmp_line = line
        flav_match = 0
        syn_match = 0
        for key in flav_list:
            if line.find(key) != -1:
                flav_match += 0
                for text in TO_REMOVE: tmp_line = tmp_line.replace(text, '&$#@!')  # replace text with unique delimiter
                tmp_list = tmp_line.split('&$#@!')  # split at delimiter
                tmp_list = list(filter(None, tmp_list))  # remove any empty elements in the list
                tmp_list[:] = [item for item in tmp_list if item != '']  # remove any empty elements in the list
                plant_flavs[plant_name].append(tmp_list)  # add the compound for the specific plant

        if flav_match == 0:
            line_mod = line.upper().replace(' ', '')
            for key in flav_synonyms:
                for flav in flav_synonyms[key]:
                    flav_mod = flav.upper().replace(' ', '')
                    if line_mod.find(flav_mod) != -1:
                        syn_match += 0
                        for text in TO_REMOVE: tmp_line = tmp_line.replace(text, '&$#@!')
                        tmp_list = tmp_line.split('&$#@!')
                        tmp_list = list(filter(None, tmp_list))
                        tmp_list[:] = [item for item in tmp_list if item != '']
                        check = re.findall(cas_id_form, tmp_list[2])  # see if the line has a CAS ID

                        # add to the output string such that the data is written in CSV format
                        if len(check) < 1: plant_flavs[plant_name].append('[Alt Name for ' + key + ']' + tmp_list[2])
                        else: plant_flavs[plant_name].append('[Alt Name for ' + key + ']' + tmp_list[3])
            if syn_match == 0:
                for key in flav_relatives:
                    for flav in flav_relatives[key]:
                        flav_mod = flav.upper().replace(' ', '')
                        if line_mod.find(flav_mod) != -1:
                            for text in TO_REMOVE: tmp_line = tmp_line.replace(text, '&$#@!')
                            tmp_list = tmp_line.split('&$#@!')
                            tmp_list = list(filter(None, tmp_list))
                            tmp_list[:] = [item for item in tmp_list if item != '']
                            check = re.findall(cas_id_form, tmp_list[2])  # see if the line has a CAS ID

                            # add to the output string such that the data is written in CSV format
                            if len(check) < 1: plant_flavs[plant_name].append('[Relative of ' + key + ']' +
                                                                              tmp_list[2])
                            else: plant_flavs[plant_name].append('[Alt Name for ' + key + ']' + tmp_list[3])
    file.close()
    # os.remove(file_name)

if __name__ == '__main__':
    main()

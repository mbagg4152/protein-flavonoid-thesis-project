from jsondata import *
from pathstrings import SEP as slash
from miscvals import *
import os
import re

new_dir = '..' + SEP + 'misc_output' + SEP + 'knapsack_dir'
ks_data_name = new_dir + SEP + 'knapsack_data.csv'
plant_flavs = {}

# sections of text to be removed from the fetched HTML file
to_remove = ['<tr>', '</tr>', '</font>', "\n", "\\n", ">", '</a>', '</td>', "target=\"_blank\">", "target=\"_blank\"",
             "<td class=\"d1\">", '<font color=#FF00BF>', '<a href=information.php?word=']

def main():
    try: os.mkdir(new_dir)
    except FileExistsError: pass

    # parsing file of plant codes and names
    for key in plant_dict_reg:
        tmp_val = plant_dict_reg[key]  # plant name
        tmp_url = KS_URL + tmp_val + "'"  # fill out URL for the plant
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
            cas_id_form = "([0-9]*-[0-9]*-[0-9]*)"  # regular expression for CAS ID
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
        for flav in flav_list:
            if line.find(flav) != -1:
                for text in to_remove: tmp_line = tmp_line.replace(text, '&$#@!')  # replace text with unique delimiter
                tmp_list = tmp_line.split('&$#@!')  # split at delimiter
                tmp_list = list(filter(None, tmp_list))  # remove any empty elements in the list
                tmp_list[:] = [item for item in tmp_list if item != '']  # remove any empty elements in the list
                plant_flavs[plant_name].append(tmp_list)  # add the compound for the specific plant
    file.close()
    os.remove(file_name)

if __name__ == '__main__':
    main()

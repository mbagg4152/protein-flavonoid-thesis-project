# import wget
from jsondata import *
from pathstrings import slash as slash
import re
import os

partial_url = "'http://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word="
cwd = os.getcwd()
new_dir = '..' + slash + 'misc_output' + slash + 'knapsack_dir'
ks_data_name = new_dir + slash + 'knapsack_data.txt'
wget_out = new_dir + slash + 'wget-out.txt'
plant_flavs = {}
to_remove = ['<tr>', '</tr>', '<a href=information.php?word=', '</a>', "target=\"_blank\">", "target=\"_blank\"",
             '</td>', "<td class=\"d1\">", '<font color=#FF00BF>', '</font>', "\n", "\\n", ">"]

def main():
    try: os.mkdir(new_dir)
    except FileExistsError: pass
    for key in plant_pairs:
        tmp_val = plant_pairs[key]
        tmp_url = partial_url + tmp_val + "'"
        tmp_file_name = new_dir + slash + key + ".csv"
        wget_download(tmp_file_name, tmp_url)
        parse_file(tmp_file_name, tmp_val)
    try:
        tmp_file = open(ks_data_name, 'x')
        tmp_file.close()
    except FileExistsError or PermissionError: print('')

    ks_file = open(ks_data_name, 'w')
    ks_str = ''
    for key in plant_flavs:
        ks_str += key + ': '
        if len(plant_flavs[key]) > 0:
            # ks_str += "\n"
            cas_id_form = "([0-9]*-[0-9]*-[0-9]*)"
            for item in plant_flavs[key]:
                check = re.findall(cas_id_form, item[2])
                if len(check) < 1: ks_str += ' ,' + item[0] + ',' + item[2] + '\n'
                else: ks_str += ' ,' + item[0] + ',' + item[3] + '\n'

        ks_str += '\n'
    ks_file.write(ks_str)
    ks_file.close()


def wget_download(name, url):
    try:
        tmp_file = open(name, 'x')
        tmp_file.close()
    except FileExistsError: pass
    # os.system("echo >" + wget_out)
    tmp_cmd = "wget -O " + name + ' ' + url
    os.system(tmp_cmd)


def parse_file(file_name, plant_name):
    file = open(file_name, 'r')
    lines = file.read().split('<tr')
    plant_flavs[plant_name] = []
    for line in lines:
        tmp_line = line
        for flav in flav_list:
            if line.find(flav) != -1:
                for text in to_remove:
                    tmp_line = tmp_line.replace(text, '&$#@!')
                tmp_list = tmp_line.split('&$#@!')
                tmp_list = list(filter(None, tmp_list))
                tmp_list[:] = [item for item in tmp_list if item != '']
                plant_flavs[plant_name].append(tmp_list)
    file.close()
    os.remove(file_name)

if __name__ == '__main__':
    main()

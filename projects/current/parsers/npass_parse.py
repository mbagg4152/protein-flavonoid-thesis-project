import re
from parsersconstants import *

plant_dict_output = ''
plant_codes = ''
comp_dict_out = ''
comp_codes_out = ''
out_list = []
kegg_names = []
kegg_upper = []
comp_names = []
comp_upper = []
plant_dict = get_json_data(plant_names_codes)
flavs = get_json_data(flav_names)
nc_codes = get_json_data(ncomp_codes)
nc_dict = get_json_data(ncomp_dict)
np_codes = get_json_data(nplant_codes)
np_dict = get_json_data(nplant_dict)
big_output = ''


def main():
    # setup_org_files()
    # setup_comp_files()
    do_output()


def do_output():
    file = open(pairs_in, 'r', encoding="utf8", errors='ignore')
    contents = file.readlines()
    file.close()

    for line in contents: proc_pairs(line)

    out_file = open(selected, 'w')
    out_file.write(big_output)
    out_file.close()


def proc_pairs(line):
    global big_output
    try:
        tmp_comp = re.findall(re_org_comp_id, line)[0].strip().strip('\t\r\n')
        tmp_plant = re.findall(re_org_name_pairs, line)[0].strip().strip('\t\r\n')
        if tmp_comp in nc_codes and tmp_plant in np_codes:
            plant_name = np_dict.get(tmp_plant)
            comp_name = nc_dict.get(tmp_comp)
            tmp_str = plant_name + ', ' + comp_name + '\n'
            if tmp_str not in big_output: big_output += tmp_str
    except IndexError:
        return
    pass


def setup_comp_files():
    global comp_names
    file = open(comp_in, 'r', encoding="utf8", errors='ignore')
    contents = file.readlines()
    file.close()
    comp_names = flavs

    for name in comp_names:
        name = "".join(name.split())
        comp_upper.append(name.upper())

    for line in contents: make_comp_data(line)

    for line in contents: make_plant_data(line)
    out_file = open(comp_dict, 'w')
    out_file.write(comp_dict_out)
    out_file.close()

    out_file = open(cmp_codes, 'w')
    out_file.write(comp_codes_out)
    out_file.close()


def make_comp_data(line):
    global comp_codes_out, comp_dict_out
    try:
        tmp_id = re.findall(re_comp_id, line)[0].strip().strip('\t\r\n')
        tmp_name = re.findall(re_comp_name, line)[0].strip().strip('\t\r\n')
        tmp_nosp = "".join(tmp_name.split())

        for comp in comp_upper:
            if comp.strip() in tmp_nosp.upper().strip() or tmp_nosp.upper().strip() in comp.strip():
                print(tmp_nosp.upper())
                tmp_str = '\"' + tmp_id + '\":\"' + tmp_name + '\",\n'

                if tmp_id not in comp_dict_out: comp_dict_out += tmp_str
                tmp_code = '\"' + tmp_id + '\",\n'
                if tmp_id not in comp_codes_out: comp_codes_out += tmp_code
    except IndexError:
        return


def setup_org_files():
    global kegg_names, kegg_upper
    file = open(org_in, 'r', encoding="utf8", errors='ignore')
    contents = file.readlines()
    file.close()
    kegg_names = list(plant_dict.values())

    for name in kegg_names:
        name = "".join(name.split())
        kegg_upper.append(name.upper())

    for line in contents: make_plant_data(line)
    out_file = open(org_out, 'w')
    out_file.write(plant_dict_output)
    out_file.close()

    out_file = open(org_codes, 'w')
    out_file.write(plant_codes)
    out_file.close()


def make_plant_data(line):
    global plant_dict_output, plant_codes
    try:
        tmp_id = re.findall(re_org_id, line)[0].strip().strip('\t\r\n')
        tmp_name = re.findall(re_org_name, line)[0].strip().strip('\t\r\n')
        tmp_nosp = "".join(tmp_name.split())

        if tmp_nosp.upper() in kegg_upper:
            print(tmp_nosp.upper())
            tmp_str = '\"' + tmp_id + '\":\"' + tmp_name + '\",'

            plant_dict_output += tmp_str
            tmp_code = '\"' + tmp_id + '\",'
            plant_codes += tmp_code
    except IndexError:
        return


if __name__ == '__main__':
    main()

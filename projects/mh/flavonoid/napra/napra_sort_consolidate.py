import numpy as np
import os
import pandas as pd

sp_dict = {}
codes = ['AGI', 'BUN', 'KXN', 'HWB', 'ERD', 'GEN', 'KMP', 'LU2', 'MYC', 'NAR', 'QUE']

full_data = []


class spec(object):
    name = ''
    family = ''
    comps = []

    def __init__(self, name, family, comps):
        self.name = name
        self.family = family
        self.comps = comps


def main():
    global codes, sp_dict, full_data
    try: os.mkdir('napra')
    except (FileExistsError, OSError) as e: pass
    rd = pd.read_excel('vrmn.xlsx', sheet_name='SPL')
    for item in codes:
        full_data = []
        tmp = rd[[item]].values.tolist()
        flat = [item for sublist in tmp for item in sublist]
        stripped = [x for x in flat if isinstance(x, str)]
        sp_dict[item] = stripped

        tmp_df = pd.read_excel('vrmn.xlsx', sheet_name=item)
        col_list = np.arange(0, tmp_df.shape[1]).tolist()
        tmp_df.columns = col_list

        for name in stripped:
            selected = tmp_df.loc[tmp_df[1] == name]
            selected = selected.reset_index()
            selected = selected.drop(columns=['index'])
            try:
                fam = selected.at[0, 0]
                idx_vals = selected.index.values.tolist()
                if (len(idx_vals) > 0) and (len(idx_vals) < 2):
                    selected = selected.iloc[0].values.tolist()[2:]
                    sel = [x for x in selected if isinstance(x, str)]
                    sel = sorted(set(sel))
                    ts = spec(name=name, family=fam, comps=sel)
                    full_data.append(ts)
                else:
                    selected = selected.drop(columns=[0, 1])
                    data = []
                    for index, row in selected.iterrows():
                        data += row.values.tolist()
                    sel = [x for x in data if isinstance(x, str)]
                    sel = sorted(set(sel))
                    ts = spec(name=name, family=fam, comps=sel)
                    full_data.append(ts)

            except KeyError: pass
        out = ''
        for fd in full_data:
            this_str = fd.family + '#' + fd.name + '#' + '#'.join(fd.comps) + '\n'
            out += this_str
        with open('napra' + os.sep + item + '.csv', "w") as text_file:
            text_file.write(out)


if __name__ == '__main__':
    main()

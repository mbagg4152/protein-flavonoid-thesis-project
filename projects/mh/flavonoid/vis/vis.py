from matplotlib.font_manager import FontProperties
from matplotlib_venn import venn2 as v2
from matplotlib_venn import venn2_circles as v2c
from matplotlib_venn import venn2_unweighted as v2u
from matplotlib_venn import venn3 as v3
from matplotlib_venn import venn3_circles as v3c
from pandas.plotting import table
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

AK = 10
AL = 21
AS = 32
K = 43
L = 54
LK = 65
S = 76
SK = 87
TP = 98
TC = 109

fname_cats = {AK: 'is_predicted_', AL: 'is_exp_known_', AS: 'is_similar_',
              K: 'only_predicted_', L: 'only_exp_known_', S: 'only_similar_',
              LK: 'known_predicted_', SK: 'similar_predicted_'}
flav_ids = ['AGI', 'BUN', 'KXN', 'HWB', 'EKXN', 'EGT', 'ERD', 'GKXN', 'GEN',
            'HCC', 'KMP', 'LU2', 'MYC', 'NAR', 'QUE']

names = {'AGI': 'Apigenin', 'BUN': 'Butein', 'KXN': 'Catechin', 'HWB': 'Cyanidin',
         'EKXN': 'Epicatechin', 'EGT': 'Epigallocatechin', 'ERD': 'Eriodictyol',
         'GKXN': 'Gallocatechin', 'GEN': 'Genistein', 'HCC': 'Isoliquiritigenin',
         'KMP': 'Kaempferol', 'LU2': 'Luteolin', 'MYC': 'Myricetin', 'NAR': 'Naringenin',
         'QUE': 'Quercetin'}

fe = fm.FontEntry(fname='nunito/NunitoSans-SemiBold.ttf', name='NunitoSans-Bold')
fm.fontManager.ttflist.insert(0, fe)
plt.rcParams['font.family'] = fe.name
ftype = '.png'

def main():
    rd = pd.read_excel('kegg-data.xlsx')
    df = rd[['Name'] + flav_ids]
    df = df.set_index('Name')
    mk_tables(df)
    # vennd(df)

def mk_tables(df):
    try: os.mkdir('pics')
    except (FileExistsError, OSError) as e: pass
    try: os.mkdir('csv')
    except (FileExistsError, OSError) as e: pass
    
    for flav in flav_ids:
        x_kegg = df.loc[df[flav] == 'K'][flav].index.values.tolist()
        x_lit = df.loc[df[flav] == 'L'][flav].index.values.tolist()
        x_lit_kegg = df.loc[df[flav] == 'LK'][flav].index.values.tolist()
        x_sim = df.loc[df[flav] == 'S'][flav].index.values.tolist()
        x_sim_kegg = df.loc[df[flav] == 'SK'][flav].index.values.tolist()
        
        data = {AK: x_kegg + x_lit_kegg + x_sim_kegg, AL: x_lit + x_lit_kegg,
                AS: x_sim_kegg + x_sim, K: x_kegg, L: x_lit, LK: x_lit_kegg,
                SK: x_sim_kegg}
        for key in data:
            if len(data[key]):
                split_list = np.array_split(np.array(name_frmt(data[key])), 4)
                tmp_df = format_df(split_list, max_lst_len(split_list), flav)
                output_table(tmp_df, get_fname(flav, key, TP))
                tmp_df.to_csv(index=False, header=False, path_or_buf=get_fname(flav, key, TC))

def get_fname(flav, cat, tab_type):
    f_ext = ''
    main_dir = ''
    if tab_type == TP:
        main_dir = 'pics'
        f_ext = '.png'
    else:
        main_dir = 'csv'
        f_ext = '.csv'
    
    try: os.mkdir(main_dir + os.sep + names[flav])
    except (FileExistsError, OSError) as e: pass
    
    return main_dir + os.sep + names[flav] + os.sep + fname_cats[cat] + flav.lower() + f_ext

def name_frmt(lst):
    lst = sorted(set(lst))
    cpy = lst.copy()
    offset = 0
    num = 1
    for idx in range(len(lst)):
        tmp = lst[idx]
        if len(lst[idx].split()) > 2:
            wrd = lst[idx].split()
            sp = ' ' * ((2 * len(str(num))) + 3)
            tmp = wrd[0] + ' ' + wrd[1]
        lst[idx] = str(num) + '. ' + tmp
        num += 1
    return lst

def max_lst_len(lst):
    if len(lst) == 5: return max(len(lst[0]), len(lst[1]), len(lst[2]), len(lst[3]), len(lst[4]))
    elif len(lst) == 4: return max(len(lst[0]), len(lst[1]), len(lst[2]), len(lst[3]))
    elif len(lst) == 3: return max(len(lst[0]), len(lst[1]), len(lst[2]))
    elif len(lst) == 2: return max(len(lst[0]), len(lst[1]))
    elif len(lst) == 1: return len(lst[0])
    else: return 1

def output_table(df, fname):
    if len(df):
        row_cnt = len(df.values.tolist())
        plt.clf()
        
        tab = plt.table(cellText=df.values, loc='center', cellLoc='left', colLoc='left',
                        edges='open')
        plt.axis('off')
        plt.axis('tight')
        tab.auto_set_font_size(False)
        vals = [0] * row_cnt
        for key, cell in tab.get_celld().items():
            if len(cell._text.get_text().split()) > 3: vals[key[0]] += 1
            cell.set_fontsize(6)
            cell.set_text_props(linespacing=1, ha='left', wrap=False)
            cell.set_height(0.03)
            cell.PAD = 0.085
            cell.set_linewidth(0)
        
        for col in range(len(df.columns)): tab.auto_set_column_width(col)
        plt.savefig(fname, bbox_inches='tight', dpi=400, format='png', pad_inches=0)

def format_df(lst, ml, comp):
    df = pd.DataFrame({'del': np.arange(ml)})
    cnt = 1
    for b in lst:
        if isinstance(b, str): b = np.array([b])
        if len(b.tolist()) == 0: b = [' ']
        else: b = b.tolist() + mk_blank_list(ml - len(b.tolist()))
        
        df['C' + str(cnt)] = b
        cnt += 1
    df = df.drop(columns=['del'])
    return df

def mk_blank_list(num):
    em = []
    for idx in range(num): em.append(' ')
    return em

def vennd(df):
    try: os.mkdir('pvenn')
    except (FileExistsError, OSError) as e: pass
    
    left = '10'
    right = '01'
    mid = '11'
    
    for comp in flav_ids:
        kvals = df.loc[df[comp] == 'kegg'][comp].index.values.tolist()
        lvals = df.loc[df[comp] == 'lit'][comp].index.values.tolist()
        bvals = df.loc[df[comp] == 'lit-kegg'][comp].index.values.tolist()
        
        if len(lvals) == 0: lvals = ['None']
        if len(bvals) == 0: bvals = ['None']
        if len(kvals) == 0: kvals = ['None']
        
        kcnt = len(kvals)
        lcnt = len(lvals)
        bcnt = len(bvals)
        
        v = v2(subsets=(lcnt, kcnt, bcnt), set_labels=('', ''))
        v.get_label_by_id(left).set_text(lcnt)
        v.get_label_by_id(mid).set_text(bcnt)
        v.get_label_by_id(right).set_text(kcnt)
        plt.annotate('EXP only', xytext=(-70, -70), ha='center', textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='gray'),
                     xy=v.get_label_by_id(left).get_position() - np.array([0, 0.05]))
        plt.annotate('KEGG only', xytext=(70, -70), ha='center', textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                     xy=v.get_label_by_id(right).get_position() - np.array([0, 0.05]),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',
                                     color='gray'))
        plt.savefig('pvenn/' + comp + '.png')

if __name__ == '__main__':
    main()

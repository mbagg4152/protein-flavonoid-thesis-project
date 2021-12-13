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

flav_ids = ['AGI', 'BUN', 'KXN', 'HWB', 'EC', 'EGT', 'ERD', 'GC', 'GEN', 'HCC', 'KMP', 'LU2', 'MYC', 'NAR', 'QUE']

names = {'AGI': 'Apigenin', 'BUN': 'Butein', 'KXN': 'Catechin', 'HWB': 'Cyanidin', 'EC': 'Epicatechin',
         'EGT': 'Epigallocatechin', 'ERD': 'Eriodictyol', 'GC': 'Gallocatechin', 'GEN': 'Genistein',
         'HCC': 'Isoliquiritigenin', 'KMP': 'Kaempferol', 'LU2': 'Luteolin', 'MYC': 'Myricetin', 'NAR': 'Naringenin',
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
    try: os.mkdir('tables')
    except (FileExistsError, OSError) as e: pass
    for flav in flav_ids:
        kd = df.loc[df[flav] == 'kegg'][flav].index.values.tolist()
        lk = df.loc[df[flav] == 'lit-kegg'][flav].index.values.tolist()
        rk = df.loc[df[flav] == 'rel-kegg'][flav].index.values.tolist()
        ld = df.loc[df[flav] == 'lit'][flav].index.values.tolist()
        rel = df.loc[df[flav] == 'rel'][flav].index.values.tolist()
        rkc = rk.copy()

        konly = name_wrap(kd)

        # for i in range(len(rk)): rk[i] = rk[i] + '**'
        # for i in range(len(rk)): rkc[i] = rkc[i] + '**'
        # for i in range(len(rel)): rel[i] = rel[i] + '**'

        lonly = name_wrap(ld + rel)
        lit = name_wrap(ld + lk)

        both = name_wrap(lk + rk)
        # for i in range(len(lk)): lk[i] = lk[i] + '*'
        kegg = name_wrap(kd + lk + rk)

        konly = np.array_split(np.array(konly), 5)
        lonly = np.array_split(np.array(lonly), 5)
        both = np.array_split(np.array(both), 5)
        rkegg = np.array_split(np.array(rk), 5)
        kegg = np.array_split(np.array(kegg), 5)
        lit = np.array_split(np.array(lit), 5)

        output_table(format_df(kegg, max_lst_len(kegg), flav, 'k'), 'tables/' + 'all_kegg_' + flav + ftype)
        output_table(format_df(konly, max_lst_len(konly), flav, 'ko'), 'tables/' + 'just_kegg_' + flav + ftype)

        output_table(format_df(lonly, max_lst_len(lonly), flav, 'lo'), 'tables/' + 'just_lit_' + flav + ftype)
        output_table(format_df(lit, max_lst_len(lit), flav, 'l'), 'tables/' + 'all_lit_' + flav + ftype)
        output_table(format_df(both, max_lst_len(both), flav, 'b'), 'tables/' + 'just_lit_kegg_' + flav + ftype)
        output_table(format_df(rkegg, max_lst_len(rkegg), flav, 'rk'), 'tables/' + 'just_rel_kegg_' + flav + ftype)


def name_wrap(lst):
    lst = sorted(set(lst))
    cpy = lst.copy()
    offset = 0
    num = 1
    for idx in range(len(lst)):
        if len(lst[idx].split()) > 2:
            wrd = lst[idx].split()
            sp = ''
            for i in range((2 * len(str(num))) + 3): sp += ' '
            # lst[idx] = str(num) + '. ' + wrd[0] + ' ' + wrd[1] + '\n' + sp + \
            #            ' '.join([x.lower() for x in wrd[2:]])
            lst[idx] = str(num) + '. ' + wrd[0] + ' ' + wrd[1]
        else: lst[idx] = str(num) + '. ' + lst[idx]
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
        tab = plt.table(cellText=df.values, loc='center', cellLoc='left', colLoc='left', edges='open')
        # tab.scale(1, 1.5)
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

        # tcells = tab.get_celld()
        # for i in range(0, len(vals)):
        #     for j in range(0, len(df.head())):
        #         tcells[(i, j)].set_fontsize(6)
        #         # if vals[i]:
        #         #     tcells[(i, j)].set_height(0.048)
        #         # else:
        #         #     tcells[(i, j)].set_height(0.048)
        #         tcells[(i, j)].set_text_props(linespacing=1, ha='left', wrap=False)
        #         tcells[(i, j)].set_height(0.03)
        #         tcells[(i, j)].PAD = 0.085
        #         tcells[(i, j)].set_linewidth(0)

        for col in range(len(df.columns)): tab.auto_set_column_width(col)
        plt.savefig(fname, bbox_inches='tight', dpi=400, format='png', pad_inches=0)


def format_df(lst, ml, comp, tag):
    df = pd.DataFrame({'del': np.arange(ml)})
    cnt = 1
    for b in lst:
        if isinstance(b, str): b = np.array([b])
        if len(b.tolist()) == 0: b = [' ']
        else: b = b.tolist() + mk_blank_list(ml - len(b.tolist()))
        df['C' + str(cnt)] = b
        cnt += 1
    df = df.drop(columns=['del'])
    # df.loc[len(df.index)] = ['', '', '', '', '']

    # if (tag != 'l') and (tag != 'lo'):
    #     df.loc[len(df.index)] = ['*  = Known to synthesize', names.get(comp), '', '', '']
    #     df.loc[len(df.index)] = ['** = Known to synthesize', 'compound(s) similar to', names.get(comp), '', '']
    # df = df[1:]
    return df


def mk_blank_list(num):
    em = []
    for idx in range(num): em.append(' ')
    return em


def vennd(df):
    try: os.mkdir('venn')
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
        plt.annotate('Experimentally\nKnown\nOnly', xy=v.get_label_by_id(left).get_position() - np.array([0, 0.05]),
                     xytext=(-70, -70), ha='center', textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='gray'))
        plt.annotate('Predicted\nOnly', xy=v.get_label_by_id(right).get_position() - np.array([0, 0.05]),
                     xytext=(70, -70), ha='center', textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5', color='gray'))
        plt.savefig('venn/' + comp + '.png')


if __name__ == '__main__':
    main()

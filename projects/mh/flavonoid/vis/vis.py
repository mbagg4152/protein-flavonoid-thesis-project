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

flav_ids = ['AGI', 'BUN', 'KXN', 'HWB', 'EKXN', 'EGT', 'ERD', 'GKXN', 'GEN', 'HCC', 'KMP', 'LU2', 'MYC', 'NAR', 'QUE']

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
    try: os.mkdir('csv-tables')
    except (FileExistsError, OSError) as e: pass
    for flav in flav_ids:
        x_kegg = df.loc[df[flav] == 'K'][flav].index.values.tolist()
        x_lit_kegg = df.loc[df[flav] == 'LK'][flav].index.values.tolist()
        x_rel_kegg = df.loc[df[flav] == 'SK'][flav].index.values.tolist()
        x_lit = df.loc[df[flav] == 'L'][flav].index.values.tolist()
        x_rel = df.loc[df[flav] == 'S'][flav].index.values.tolist()

        all_kegg = np.array_split(np.array(name_wrap(x_kegg + x_lit_kegg + x_rel_kegg)), 4)
        all_lit = np.array_split(np.array(name_wrap(x_lit + x_lit_kegg)), 4)
        all_rel = np.array_split(np.array(name_wrap(x_rel_kegg + x_rel)), 4)
        excl_kegg = np.array_split(np.array(name_wrap(x_kegg)), 4)
        excl_lit_kegg = np.array_split(np.array(name_wrap(x_lit_kegg)), 4)
        excl_lit = np.array_split(np.array(name_wrap(x_lit)), 4)
        excl_sim_kegg = np.array_split(np.array(x_rel_kegg), 4)

        all_kegg_df = format_df(all_kegg, max_lst_len(all_kegg), flav)
        all_lit_df = format_df(all_lit, max_lst_len(all_lit), flav)
        all_sim_df = format_df(all_rel, max_lst_len(all_rel), flav)
        excl_kegg_df = format_df(excl_kegg, max_lst_len(excl_kegg), flav)
        excl_lit_df = format_df(excl_lit, max_lst_len(excl_lit), flav)
        excl_lit_kegg_df = format_df(excl_lit_kegg, max_lst_len(excl_lit_kegg), flav)
        excl_sim_kegg_df = format_df(excl_sim_kegg, max_lst_len(excl_sim_kegg), flav)

        output_table(all_kegg_df, 'tables/' + 'all_kegg_' + flav + ftype)  # all predicted entries
        output_table(all_lit_df, 'tables/' + 'all_lit_' + flav + ftype)  # all entries with exact lit matches
        output_table(all_sim_df, 'tables/' + 'all_sim_' + flav + ftype)  # all entries with sim lit matches (no exact)
        output_table(excl_kegg_df, 'tables/' + 'just_kegg_' + flav + ftype)  # predicted with no matches
        output_table(excl_lit_df, 'tables/' + 'just_lit_' + flav + ftype)  # exact match but no predictions
        output_table(excl_lit_kegg_df, 'tables/' + 'just_lit_kegg_' + flav + ftype)  # exact match AND predicted
        output_table(excl_sim_kegg_df, 'tables/' + 'just_sim_kegg_' + flav + ftype)  # sim match AND predicted

        all_kegg_df.to_csv(index=False, header=False, path_or_buf='csv-tables/all_kegg_' + flav + '.csv')
        all_lit_df.to_csv(index=False, header=False, path_or_buf='csv-tables/all_lit_' + flav + '.csv')
        all_sim_df.to_csv(index=False, header=False, path_or_buf='csv-tables/all_sim_' + flav + '.csv')
        excl_kegg_df.to_csv(index=False, header=False, path_or_buf='csv-tables/just_kegg_' + flav + '.csv')
        excl_lit_df.to_csv(index=False, header=False, path_or_buf='csv-tables/just_lit_' + flav + '.csv')
        excl_lit_kegg_df.to_csv(index=False, header=False, path_or_buf='csv-tables/just_lit_kegg_' + flav + '.csv')
        excl_sim_kegg_df.to_csv(index=False, header=False, path_or_buf='csv-tables/just_sim_kegg_' + flav + '.csv')


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
            if num < 10: lst[idx] = str(num) + '. ' + wrd[0] + ' ' + wrd[1]
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


def tweak_for_csv(df):
    col = df.columns
    for i in range(0, df.shape[0]):
        for j in range(0, len(col)):
            sp = (len(str((j + 2)))) * ' '
            try: df.at[i, col[j]] = df.at[i, col[j]].replace('. ', ".").replace(' ', '\n' + sp).replace('.', '. ')
            except KeyError:
                print(df)
                print(str(i) + str(j))
    return df


def output_table(df, fname):
    if len(df):
        row_cnt = len(df.values.tolist())
        plt.clf()
        tab = plt.table(cellText=df.values, loc='center', cellLoc='left', colLoc='left', edges='open')
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
        plt.annotate('Experimentally\nKnown\nOnly', xy=v.get_label_by_id(left).get_position() - np.array([0, 0.05]),
                     xytext=(-70, -70), ha='center', textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='gray'))
        plt.annotate('Predicted\nOnly', xy=v.get_label_by_id(right).get_position() - np.array([0, 0.05]),
                     xytext=(70, -70), ha='center', textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5', color='gray'))
        plt.savefig('pvenn/' + comp + '.png')


if __name__ == '__main__':
    main()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2 as v2
from matplotlib_venn import venn2_circles as v2c
from matplotlib_venn import venn2_unweighted as v2u
from matplotlib_venn import venn3 as v3
from matplotlib_venn import venn3_circles as v3c
import os
from pandas.plotting import table
from matplotlib.font_manager import FontProperties

comp_ids = ['AGI', 'BUN', 'KXN', 'HWB', 'EC', 'EGT', 'ERD', 'GC', 'GEN', 'HCC', 'KMP', 'LU2', 'MYC', 'NAR', 'QUE']
# comp_ids = ['AGI', 'BUN', 'KXN']
names = {'AGI': 'Apigenin', 'BUN': 'Butein', 'KXN': 'Catechin', 'HWB': 'Cyanidin', 'EC': 'Epicatechin',
         'EGT': 'Epigallocatechin', 'ERD': 'Eriodictyol', 'GC': 'Gallocatechin', 'GEN': 'Genistein',
         'HCC': 'Isoliquiritigenin', 'KMP': 'Kaempferol', 'LU2': 'Luteolin', 'MYC': 'Myricetin', 'NAR': 'Naringenin',
         'QUE': 'Quercetin'}
kegglbl = 'Predicted to be\ncapable of making\n'
litlbl = 'Known to make\n'
lklbl = 'Predicted, known to\nmake '
rellbl = 'Known to make\n'
rklbl = 'Predicted, known\nto make '


def main():
    rd = pd.read_excel('kegg-data.xlsx')
    df = rd[['Name'] + comp_ids]
    df = df.set_index('Name')
    tab3(df)


def populate(lst, ml, comp, tag):
    df = pd.DataFrame({'del': np.arange(ml)})
    cnt = 1
    for b in lst:
        if isinstance(b, str): b = np.array([b])
        if len(b.tolist()) == 0: b = [' ']
        else: b = b.tolist() + gen_empty(ml - len(b.tolist()))
        df[comp.lower() + '-' + tag + str(cnt) + '\n\n\n'] = b
        cnt += 1
    df = df.drop(columns=['del'])
    df.loc[len(df.index)] = ['', '', '', '', '']

    if (tag != 'l') and (tag != 'lo'):
        df.loc[len(df.index)] = ['*  = Known to synthesize', names.get(comp), '', '', '']
        df.loc[len(df.index)] = ['** = Known to synthesize', 'compound(s) similar to', names.get(comp), '', '']

    return df


def tab3(df):
    try: os.mkdir('tab3')
    except (FileExistsError, OSError) as e: pass
    for comp in comp_ids:
        kd = df.loc[df[comp] == 'kegg'][comp].index.values.tolist()
        lk = df.loc[df[comp] == 'lit-kegg'][comp].index.values.tolist()
        rk = df.loc[df[comp] == 'rel-kegg'][comp].index.values.tolist()
        ld = df.loc[df[comp] == 'lit'][comp].index.values.tolist()
        rel = df.loc[df[comp] == 'rel'][comp].index.values.tolist()
        rkc = rk.copy()

        ko = mod_list(kd)

        for i in range(len(rk)): rk[i] = rk[i] + '**'
        for i in range(len(rk)): rkc[i] = rkc[i] + '**'
        for i in range(len(rel)): rel[i] = rel[i] + '**'

        lo = mod_list(ld + rel)
        lit = mod_list(ld + lk)

        both = mod_list(lk + rk)
        for i in range(len(lk)): lk[i] = lk[i] + '*'
        kegg = mod_list(kd + lk + rk)

        ko = np.array_split(np.array(ko), 5)
        lo = np.array_split(np.array(lo), 5)
        both = np.array_split(np.array(both), 5)
        kegg = np.array_split(np.array(kegg), 5)
        lit = np.array_split(np.array(lit), 5)

        # plot_stuff(populate(ko, gml(ko), comp, 'ko'), 'tab3/' + comp.lower() + '-onlykegg' + '.png')
        plot_stuff(populate(kegg, gml(kegg), comp, 'k'), 'tab3/' + comp.lower() + '-kegg' + '.png')
        # plot_stuff(populate(lo, gml(lo), comp, 'lo'), 'tab3/' + comp.lower() + '-onlylit' + '.png')
        plot_stuff(populate(lit, gml(lit), comp, 'l'), 'tab3/' + comp.lower() + '-lit' + '.png')
        # plot_stuff(populate(both, gml(both), comp, 'b'), 'tab3/' + comp.lower() + '-both' + '.png')


def gml(lst):
    if len(lst) == 5: return max(len(lst[0]), len(lst[1]), len(lst[2]), len(lst[3]), len(lst[4]))
    elif len(lst) == 4: return max(len(lst[0]), len(lst[1]), len(lst[2]), len(lst[3]))
    elif len(lst) == 3: return max(len(lst[0]), len(lst[1]), len(lst[2]))
    elif len(lst) == 2: return max(len(lst[0]), len(lst[1]))
    elif len(lst) == 1: return len(lst[0])
    else: return 1


def plot_stuff(df, fname):
    if len(df):
        plt.figure(figsize=(6, 6))
        tab = plt.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='left',
                        colLoc='left', edges='open')
        plt.axis('off')
        tab.auto_set_font_size(False)
        for key, cell in tab.get_celld().items():
            if key[0] != 0 or key[1] != -1:
                cell.set_fontsize(5)
                cell.PAD = 0.03
                cell.set_text_props(linespacing=1)
                cell.set_height(0.018)
        for col in range(len(df.columns)): tab.auto_set_column_width(col)
        plt.savefig(fname, bbox_inches='tight', dpi=300)


def tab2(df):
    try: os.mkdir('tab2')
    except (FileExistsError, OSError) as e: pass
    for comp in comp_ids:
        kegg = df.loc[df[comp] == 'kegg'][comp].index.values.tolist()
        lit = df.loc[df[comp] == 'lit'][comp].index.values.tolist()
        lk = df.loc[df[comp] == 'lit-kegg'][comp].index.values.tolist()
        rel = df.loc[df[comp] == 'rel'][comp].index.values.tolist()
        rk = df.loc[df[comp] == 'rel-kegg'][comp].index.values.tolist()

        kegg = mod_list(kegg)
        lit = mod_list(lit)
        lk = mod_list(lk)
        rel = mod_list(rel)
        rk = mod_list(rk)
        name = names.get(comp)
        kegg_lbl = kegglbl + name
        lit_lbl = litlbl + name
        lk_lbl = lklbl + name
        rel_lbl = rellbl + name + '\nrelative(s)'
        rk_lbl = rklbl + name + '\nrelative(s)'
        maxlen = max(len(kegg), len(lit), len(lk), len(rel), len(rk))
        comp_df = pd.DataFrame({'del': np.arange(maxlen)})
        if len(kegg) > 0:
            comp_df = mod_df(kegg, comp_df, kegg_lbl, maxlen)
        if len(lit) > 0:
            comp_df = mod_df(lit, comp_df, lit_lbl, maxlen)
        if len(lk) > 0: comp_df = mod_df(lk, comp_df, lk_lbl, maxlen)
        if len(rel) > 0: comp_df = mod_df(rel, comp_df, rel_lbl, maxlen)
        if len(rk) > 0: comp_df = mod_df(rk, comp_df, rk_lbl, maxlen)
        comp_df = comp_df.drop(columns=['del'])

        plt.figure(figsize=(12, 12))
        tab = plt.table(cellText=comp_df.values, colLabels=comp_df.columns, loc='center', cellLoc='left', colLoc='left',
                        edges='open')

        # for cell in tab.get_celld().values(): cell.visible_edges = ''

        plt.axis('off')
        tab.auto_set_font_size(False)
        for key, cell in tab.get_celld().items():
            if key[0] == 0 or key[1] == -1:
                cell.set_height(0.04)
                cell.set_fontsize(12)
                cell.set_text_props(fontproperties=FontProperties(weight='bold'))
                cell.set_text_props(va='bottom')
                cell.set_text_props(linespacing=1.5)

            if key[0] != 0 or key[1] != -1:
                cell.set_fontsize(5)
                cell.PAD = 0.03
                cell.set_text_props(linespacing=1)
                cell.set_height(0.01)
        for col in range(len(comp_df.columns)):
            tab.auto_set_column_width(col)

        plt.savefig('tab2/all-' + comp.lower() + '.png', bbox_inches='tight', dpi=300)


def mod_df(lst, df, lbl, maxlen):
    if len(lst) > 70:
        chunked = list(chunk(lst, len(lst) // 3))
        df[lbl + '\n(1/3)'] = chunked[0] + gen_empty(maxlen - len(chunked[0]))
        df[lbl + '\n(2/3)'] = chunked[1] + gen_empty(maxlen - len(chunked[1]))
        df[lbl + '\n(3/3)'] = chunked[2] + gen_empty(maxlen - len(chunked[2]))
    elif 70 > len(lst) > 29:
        chunk1 = lst[:len(lst) // 2]
        chunk2 = lst[len(lst) // 2:]
        df[lbl + '\n(1/2)'] = chunk1 + gen_empty(maxlen - len(chunk1))
        df[lbl + '\n(2/2)'] = chunk2 + gen_empty(maxlen - len(chunk2))
    else:
        df[lbl] = lst + gen_empty(maxlen - len(lst))
    return df


def tab1(df):
    try: os.mkdir('tab1')
    except (FileExistsError, OSError) as e: pass

    for comp in comp_ids:
        kegg = df.loc[df[comp] == 'kegg'][comp].index.values.tolist()
        lit = df.loc[df[comp] == 'lit'][comp].index.values.tolist()
        lk = df.loc[df[comp] == 'lit-kegg'][comp].index.values.tolist()
        rel = df.loc[df[comp] == 'rel'][comp].index.values.tolist()
        rk = df.loc[df[comp] == 'rel-kegg'][comp].index.values.tolist()
        stuff = []
        for item in range(len(kegg)):
            if item == 0: stuff.append('KEGG')
            else: stuff.append(' ')
        for item in range(len(lit)):
            if item == 0: stuff.append('Literature')
            else: stuff.append(' ')
        for item in range(len(lk)):
            if item == 0: stuff.append('KEGG & Literature')
            else: stuff.append(' ')
        for item in range(len(rel)):
            if item == 0: stuff.append('Related Compound')
            else: stuff.append(' ')
        for item in range(len(rk)):
            if item == 0: stuff.append('KEGG & Related Compound')
            else: stuff.append(' ')

        data = list(chunk(kegg + lit + lk + rel + rk, 1))
        plt.figure(figsize=(12, 12))
        tab = plt.table(cellText=data, rowLabels=stuff, loc='center', cellLoc='left')
        tab.auto_set_column_width(col=list(range(len(data))))
        for cell in tab.get_celld().values(): cell.visible_edges = ''
        plt.axis('off')
        plt.savefig('tab1/' + comp + '.png', bbox_inches='tight')


def vennd(df):
    try: os.mkdir('venn')
    except (FileExistsError, OSError) as e: pass

    left = '10'
    right = '01'
    mid = '11'

    for comp in comp_ids:
        kvals = df.loc[df[comp] == 'kegg'][comp].index.values.tolist()
        lvals = df.loc[df[comp] == 'lit'][comp].index.values.tolist()
        bvals = df.loc[df[comp] == 'lit-kegg'][comp].index.values.tolist()

        if len(lvals) == 0: lvals = ['None']
        if len(bvals) == 0: bvals = ['None']
        if len(kvals) == 0: kvals = ['None']

        kcnt = len(kvals)
        lcnt = len(lvals)
        bcnt = len(bvals)

        if lcnt == 1: lcnt = 4
        if bcnt == 1: bcnt = 4
        if kcnt == 1: kcnt = 4

        klbl = ''
        blbl = ''
        llbl = ''

        if (kcnt >= (3 * lcnt)) and (kcnt >= (3 * bcnt)): klbl = ',\n'.join(chunk_str(kvals, 3))
        elif (kcnt >= (2 * lcnt)) and (kcnt >= (2 * bcnt)): klbl = ',\n'.join(chunk_str(kvals, 2))
        else: klbl = ',\n'.join(kvals)

        if (bcnt >= (3 * lcnt)) and (bcnt >= (3 * kcnt)): blbl = ',\n'.join(chunk_str(bvals, 3))
        elif (bcnt >= (2 * lcnt)) and (bcnt >= (2 * kcnt)): blbl = ',\n'.join(chunk_str(bvals, 2))
        else: blbl = ',\n'.join(bvals)

        if (lcnt >= (3 * bcnt)) and (lcnt >= (3 * kcnt)): llbl = ',\n'.join(chunk_str(lvals, 3))
        elif (lcnt >= (2 * bcnt)) and (lcnt >= (2 * kcnt)): llbl = ',\n'.join(chunk_str(lvals, 2))
        else: llbl = ',\n'.join(lvals)

        plt.figure(figsize=(12, 12))
        v = v2(subsets=(lcnt, kcnt, bcnt), set_labels=('Literature', 'Kegg'))
        v.get_label_by_id(left).set_text(llbl)
        v.get_label_by_id(mid).set_text(blbl)
        v.get_label_by_id(right).set_text(klbl)
        plt.savefig('venn/' + comp + '.png')


def mod_list(lst):
    lst = sorted(set(lst))
    cpy = lst.copy()
    offset = 0
    for idx in range(len(lst)):
        if len(lst[idx].split()) > 2:
            sp = lst[idx].split()
            mini_lst = []
            mini_lst = [' '.join([sp[0], sp[1]]), '  ' + ' '.join(sp[2:])]
            cpy = (cpy[:idx + offset]) + mini_lst + (cpy[idx + 1 + offset:])
            offset += 1
    return cpy


def gen_empty(num):
    em = []
    for idx in range(num): em.append(' ')
    return em


def chunk_str(ll, n):
    for i in range(0, len(ll), n):
        yield ', '.join(ll[i:i + n])


def chunk(ll, n):
    for i in range(0, len(ll), n):
        yield ll[i:i + n]


if __name__ == '__main__':
    main()

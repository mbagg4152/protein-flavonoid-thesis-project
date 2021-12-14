import pandas as pd
import pyupset as pu
import numpy as np
import matplotlib
import os

dbn = ['USDA', 'KNApSAcK', 'NPASS', 'Microbiome', 'IMPPAT', 'Dr. Duke', 'NAPRALERT']
raw_df = pd.read_excel('db.xlsx', index_col=0)


def main():
    try: os.mkdir('db')
    except (FileExistsError, OSError) as e: pass
    df = pd.DataFrame({dbn[0]: gcv(0), dbn[1]: gcv(1), dbn[2]: gcv(2), dbn[3]: gcv(3),
                       dbn[4]: gcv(4), dbn[5]: gcv(5), dbn[6]: gcv(6),
                       'count': [1] * len(raw_df.index)}, index=raw_df.index)

    grouped = df.groupby(dbn).count().sort_values('count')
    grouped.to_csv(path_or_buf='db/grouped.csv')


def gcv(idx): return raw_df[dbn[idx]].values


if __name__ == '__main__':
    main()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_excel('kegg-data.xlsx')
sel_df = df[
    ['Name', 'AGI', 'BUN', 'KXN', 'HWB', 'EC', 'EGT', 'ERD', 'GC', 'GEN', 'HCC', 'KMP', 'LU2', 'MYC', 'NAR', 'QUE']]
print(sel_df)

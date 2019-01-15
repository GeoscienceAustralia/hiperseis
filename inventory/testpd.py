import numpy as np
import pandas as pd

cols = ['begin', 'end', 'Z']


def valid_rows(df):
    '''Get index of valid rows matching first row'''
    key = df[cols]
    mask = (df[cols] == key)
    print("Mask:\n{0}".format(mask))
    rmask = np.all(mask, axis=1)
    print("Rmask:\n{0}".format(rmask))
    idx = df.index[rmask]
    print("Not empty: {0}".format(not idx.empty))
    return idx

df0 = pd.DataFrame({'begin': ['2012-12-13 21:30:00', ''], 'end': ['', '']})
# print(df0)
# print(df0.dtypes)
df0['begin'] = pd.to_datetime(df0['begin'])
df0['end'] = pd.to_datetime(df0['end'])
df0['Z'] = 'MHZ'
# print(df0)
# print(df0.dtypes)

df1 = pd.DataFrame({'begin': ['2012-12-13 21:30:00', ''], 'end': ['2015-10-16 18:00:00', '']})
df1['begin'] = pd.to_datetime(df1['begin'])
df1['end'] = pd.to_datetime(df1['end'])
df1['Z'] = 'MHZ'

print(valid_rows(df0))
print(valid_rows(df1))

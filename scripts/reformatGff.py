import pandas as pd
import sys
df = pd.read_csv(sys.argv[1],header=None)
df.columns = ['column_1']
rows = df.loc[0:1].copy()
df2 = df.loc[2:].copy()
df2['ID'] = range(1,len(df2)+1)
df2['ID'] = 'ID=' + df2['ID'].astype(str)
df2['column_1'] = df2['column_1'].str.cat(df2['ID'],sep=";")
df2 = rows.append(df2)
df2 = df2.drop( columns='ID')
df2.to_csv(sys.argv[2], index=False,header=False)
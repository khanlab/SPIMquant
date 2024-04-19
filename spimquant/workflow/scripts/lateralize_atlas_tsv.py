import pandas as pd

df_left = pd.read_csv(snakemake.input.tsv,sep='\t')
df_right = pd.read_csv(snakemake.input.tsv,sep='\t')

df_left['name'] = 'left ' + df_left['name']
df_right['name'] = 'right ' + df_right['name']

df_left['abbreviation'] = df_left['abbreviation'] + '-L'
df_right['abbreviation'] = df_right['abbreviation'] + '-R'

df_right['index' ] = df_right['index'] + 10000

pd.concat((df_left,df_right)).to_csv(snakemake.output.tsv,sep='\t',index=False)

import pandas as pd
import glob,os

sample_file=str(snakemake.input)

df_g1=pd.read_csv(sample_file, sep='\t', usecols = ['#gene1'])
df_g2=pd.read_csv(sample_file, sep='\t', usecols = ['gene2'])

df_g1=df_g1.rename(columns={'#gene1': 'gene'})
df_g2=df_g2.rename(columns={'gene2': 'gene'})

df=df_g1.append(df_g2)

df['oncoprint_type']='FUS_mate'

df=df.drop_duplicates(subset=['gene'], keep='last')

with open(snakemake.output[0], 'w') as output_f:
    print(df.to_csv(sep="\t", index=False), file=output_f)
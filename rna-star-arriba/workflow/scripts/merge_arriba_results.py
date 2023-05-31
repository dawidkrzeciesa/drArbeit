import pandas as pd
import glob,os

input_files=str(snakemake.input)

input_files= input_files.split()

df = pd.concat([pd.read_csv(f,sep="\t").assign(SampleID=os.path.basename(f).split('.')[0]) for f in input_files])
col_name="SampleID"
first_col = df.pop(col_name)
df.insert(0, col_name, first_col)

with open(snakemake.output[0], 'w') as output_f:
    print(df.to_csv(sep="\t", index=False), file=output_f)

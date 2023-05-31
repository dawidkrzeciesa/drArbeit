import pandas as pd
import glob,os

input_files=str(snakemake.input)

input_files= input_files.split()

df = pd.concat([pd.read_csv(f,sep="\t").assign(SampleID=os.path.basename(f).split('.')[0]) for f in input_files])


df=df.pivot(index="SampleID",columns="gene")

# "replace header" after pivot
df.columns=[col[1] for col in df.columns]


with open(snakemake.output[0], 'w') as output_f:
    print(df.to_csv(sep="\t", index=True), file=output_f)
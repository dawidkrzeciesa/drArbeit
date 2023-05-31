import pandas as pd
import glob,os

input_files=str(snakemake.input)

input_files= input_files.split()

df = pd.concat([pd.read_csv(f,sep="\t").assign(SampleID=os.path.basename(f).split('.')[0]) for f in input_files])

df["oncoprint_type"]="FUS"

df=df[["SampleID","#gene1","gene2","oncoprint_type"]]

df["fusions_temp"]=df["SampleID"]+'_'+df['#gene1']+'_'+df['gene2']

df=df.drop_duplicates(subset=['fusions_temp'], keep='last')

df["fusions"]=df['#gene1']+'__'+df['gene2']

df=df[["SampleID","fusions","oncoprint_type"]]

df=df.pivot(index="SampleID",columns="fusions")

# "replace header" after pivot
df.columns=[col[1] for col in df.columns]


with open(snakemake.output[0], 'w') as output_f:
    print(df.to_csv(sep="\t", index=True), file=output_f)
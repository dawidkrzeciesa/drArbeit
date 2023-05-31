import pandas as pd
import glob,os

input_files=str(snakemake.input)

input_files= input_files.split()

counts = [pd.read_csv(f, index_col=0, usecols=[0,3], sep="\t", header=None) for f in input_files]

matrix = pd.concat(counts, axis=1)

SampleID=[os.path.dirname(f).split('/')[-1] for f in input_files]

matrix.columns = [SampleID]
matrix = matrix.reset_index(level=0)
matrix.rename(columns = {0:"gene_id"}, inplace = True)

with open(snakemake.output[0], 'w') as output_f:
    print(matrix.to_csv(sep="\t", index=False), file=output_f)
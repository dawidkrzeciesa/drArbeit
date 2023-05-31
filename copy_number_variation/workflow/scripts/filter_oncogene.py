import pandas as pd
import glob,os


####### OncoKB #######

inport_columns=["Hugo Symbol", "Is Oncogene", "Is Tumor Suppressor Gene"]

OncoKB = pd.read_csv(snakemake.input["oncokb"], sep="\t", usecols=inport_columns)

OncoKB = OncoKB.rename(columns={"Hugo Symbol": "gene", "Is Oncogene": "oncogene", "Is Tumor Suppressor Gene": "tsg"})

OncoKB = OncoKB[OncoKB.oncogene != "No"]

oncogene_lst = OncoKB["gene"].tolist()

####### filter cns #######
cns = pd.read_csv(snakemake.input["cns"], sep="\t")

cns["gene"]=cns.gene.str.split(",")
cns = cns.explode("gene")

cns_oncogene = cns[cns["gene"].isin(oncogene_lst)]

cns_oncogene.to_csv(snakemake.output["cns_oncogene"], sep="\t", index=False)

####### filter cnr #######

cnr = pd.read_csv(snakemake.input["cnr"], sep="\t", low_memory=False)

cnr_oncogene = cnr[cnr["gene"].isin(oncogene_lst)]

cnr_oncogene.to_csv(snakemake.output["cnr_oncogene"], sep="\t", index=False)

import os
import json
import pandas as pd

results=[]

for file in os.listdir():
    if file.endswith('result.tsv'):
        results.append(file)

ids=[]
alleles=[]
for file in results:
    with open(file,'r') as file:
        sample = file.name.split('.')[0][:-7]
        ids.append(sample)
        genotype=file.read()
        haplotype=genotype.split('\t')[-8:-2]
        alleles.append(haplotype)

hla = dict(zip(ids, alleles))
df=pd.DataFrame.from_dict(hla,orient='index',columns=['A1','A2','B1','B2','C1','C2'])
df.to_csv('hla.tsv',sep="\t")

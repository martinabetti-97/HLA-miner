import pandas as pd

#-----------------------------importing data-------------------------------

hla = pd.read_csv("hla.tsv",sep='\t',header=0,index_col=0)
info = pd.read_csv("info.tsv",sep='\t',header=0,index_col=0)
st = pd.read_excel("supertypes.xls",sep=';',header=4,index_col=1)

#---------------------------hla&info modification--------------------------
print("Setting simplified sample IDs")
#exlcuding na values ids
info = info.dropna(axis=0)
merged = hla.merge(info,how="inner",right_index=True,left_index=True)
samples = merged.index

#new ids creation
c=0
i=[]
for p in samples:
    c+=1
    p = f'p{c}'
    i.append(p)

#applying new ids to df
corr=pd.DataFrame([samples,i])
hla = hla.loc[samples, :]
info = info.loc[samples, :]
hla.index=i
info.index=i
merged.index=i
info.index.name='samples'

#exporting data
hla.to_csv('hla.tsv',sep='\t')
info.to_csv('info.tsv',sep='\t')
corr.to_csv('ids.tsv',sep='\t')

#--------------------------------supertypes--------------------------------
print("Assigning supertypes")
#ST assignment
st['Supertype'] = st['Supertype'].str[:3]
supert=st['Supertype'].tolist()
alleles=st.index.tolist()
hlaST = hla.replace({':':''}, regex=True) 
d = dict(zip(alleles, supert))
for allele in d.keys():
        hlaST=hlaST.replace(d)
hlaST=hlaST.replace({'Unc':'NC'}, regex=True)
hlaST=hlaST.replace({' ':'/'},regex=True)

#exlcuding non common ids
mergedST = hlaST.merge(info,how="inner",right_index=True,left_index=True)

#exporting data
hlaST.to_csv(f'hlaST.tsv', sep="\t")

#-------------------------division by respose group------------------------
print("Diving samples into subgroups")
NR_hla = merged.loc[merged['group'] =="NR"]
R_hla = merged.loc[merged['group'] == "R"]
NR_hlaST = mergedST.loc[mergedST['group'] == "NR"]
R_hlaST = mergedST.loc[mergedST['group']== "R"]

#exporting data
R_hla.iloc[:,:6].to_csv('hla_R.tsv',sep='\t')
NR_hla.iloc[:,:6].to_csv('hla_NR.tsv',sep='\t')
R_hlaST.iloc[:,:6].to_csv('hlaST_R.tsv',sep='\t')
NR_hlaST.iloc[:,:6].to_csv('hlaST_NR.tsv',sep='\t')

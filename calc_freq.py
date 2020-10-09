import pandas as pd
import sys
import os
import json
import argparse
import seaborn as sns
from matplotlib import pyplot as plt

#optional ethnicity analysis
e=sys.argv[1]

#----------------------relative frequencies calculation--------------------
def calc_freq(gene,df):
    col1 = df[f'{gene}1'].values.tolist()
    col2 = df[f'{gene}2'].values.tolist() 
    a = col1 + col2
    n = len(a)
    freq = {x:a.count(x)/n for x in a}
    return freq

def export(f,group):
    df = pd.read_csv(file, sep='\t', header = 0,index_col=0)
    d1 =  {}
    if df.empty: #no individual in the dataset belonging to this group
        with open(f'freq_{file}', 'w') as fp: 
            fp.write('alleles\tfrequency \n Nan\tNan')
    else:
        for gene in ['A','B','C']:
            d2 = calc_freq(gene,df)
            d1.update(d2)
        freq_df = pd.DataFrame.from_dict(data=d1, orient='index')
    freq_df.columns=['frequency']
    freq_df.to_csv(f'freq_{group}.tsv',sep="\t")
    return 1

#--------------------applying frequency calculation function---------------

subg=[]
path = os.getcwd()

for file in os.listdir(path):
    if file.startswith('hla'):
        subg.append(file)

for f in subg:
    with open(f,'r') as file:
        group = file.name.split('.')[0]
        export(file,group)

#----------------------------merging results-------------------------------

#groups to merge
et=['asian','arab','austronasian','black','caucasian']
re=['NR','R']

#calculation of difference in frequency between response groups
def difference(df):
    for a in df:
        df['diff'] = df['R'].sub(df['NR'], axis = 0)
        df['diff'] = df['diff'].abs()
    return df

#merge intermediate frequency files on allele names
#if an allele is not present in a file, value 0 is given
def merge(in_filename,groups,out_filename):
    df=pd.read_csv(in_filename,sep='\t', header= 0)
    df.columns=['alleles','global']
    for e in groups: 
        freq=pd.read_csv(f'{in_filename[:-4]}_{e}.tsv',sep='\t', header=1)
        freq.columns=['alleles','frequency']
        freq.dropna(axis=0)
        alleles=[x.strip(' ') for x in freq['alleles']] #bug 
        freq['alleles']=alleles
        freq=freq.set_index('alleles')
        freq.columns=[e]
        df=df.join(freq,how='outer',on='alleles')
        df.fillna(0, inplace=True) #freq is set to zero when a groups does not have any occurence for a specific allele
    if groups == re:
        df = difference(df)
    df=df.set_index('alleles',drop=True)
    if out_filename == "resp_st":
        for st in df.index:
            if '*' in str(st) or 'N' in str(st): #C genes and ND rows index
                df=df.drop(st,axis=0)
    df.to_csv(f'{out_filename}_freq.tsv', sep='\t')
    return df

#two frequency tables are generted, one for response group by allele and one for response gorup by supertype 
merge("freq_hla.tsv",re,"resp_a")
merge("freq_hlaST.tsv",re,"resp_st")

#frequency table for ethnic groups (optional)
if e == 'on':
    merge("freq_hla.tsv",et,"ethnicity")
    
#-------------------------------presence/absence---------------------------
def counts(df,group):
    d={}
    for gene in ['A','B','C']:
        col1 = df[f'{gene}1'].values.tolist()
        col2 = df[f'{gene}2'].values.tolist() 
        d1 = dict.fromkeys(col1+col2,0)
        for i in range(len(col1)):
            d1[col1[i]]=d1[col1[i]]+1
            if col1[i] != col2[i]:
                d1[col2[i]]=d1[col2[i]]+1   
        d.update(d1)    
    df = pd.DataFrame.from_dict(data=d, orient='index')
    df.columns=[group]
    return df

def export(out,input):
    R_all = pd.read_csv(f"{input}_R.tsv", sep='\t', header = 0,index_col=0)
    NR_all = pd.read_csv(f"{input}_NR.tsv", sep='\t', header = 0,index_col=0)
    R=counts(R_all,'R')
    NR=counts(NR_all,'NR')
    df=R.join(NR,how='outer')
    df.fillna(0, inplace=True) 
    df.to_csv(f'{out}.tsv', sep='\t')
    return df

export("count_st","hlaST")
export("count_a","hla")

#------------------------------homozygosity--------------------------------

# check homozygosity to save counts and alleles of interest and calculate frequencies
def homozygosity(df,sg):
    alleles=[]
    fdic={}
    cdic={}
    for gene in ['A','B','C']:
        count=0
        col1=list(df[f'{gene}1'])
        col2=list(df[f'{gene}2'])
        if len(col1) == 0:
            pass
        else:
            for i in range(len(col1)):
                if col1[i] == col2[i]:
                    count+=1
                    alleles.append(col1[i]) #list of alleles involved in a homozygosity event
            freq = count/len(col1)
        fdic[gene]=freq # dictionary frequency homozigous
        cdic[gene]=[count,len(col1)-count] # dictionary count of homozigous vs heterozygous
    return fdic,cdic,alleles

#extra info files
def count_alleles(all_dict):
    NR = all_dict['NR']
    R = all_dict['R']
    d={}
    for e in list(dict.fromkeys(NR+R)):
        if len(NR) != 0:
            NR_c = NR.count(e)/len(NR) #normalized counts
        else:
            NR_c = 0
        if len(R) != 0:
            R_c = R.count(e)/len(R) #normalized counts
        else:
            R_c = 0
        d[e]=[NR_c, R_c] 
    alleles_df= pd.DataFrame.from_dict(d,orient='index', columns=['NR','R'])
    return alleles_df

#merge the frequency obtained for different groups in a single table 
def export_objects():
    counts=[]
    all_dict={}
    freq_df = pd.DataFrame(columns=['NR','R'],index=['A','B','C'])
    for sg in ['NR','R']:
        df = pd.read_csv(f"hla_{sg}.tsv", sep='\t', header = 0)
        freq, count, alleles = homozygosity(df,sg)
        freq_df[sg]=freq.values()
        counts.append(count)
        all_dict[sg]=alleles
    for gene in ['A','B','C']:
        g = pd.DataFrame([counts[0][gene],counts[1][gene]],columns=['homo','hetero'],index=['NR','R'])
        g.to_csv(f'{gene}_homoCounts.tsv', sep="\t")
    alleles_df=count_alleles(all_dict)
    freq_df.to_csv('homo_freq.tsv', sep="\t")
    alleles_df.to_csv('homo_alleleCounts.tsv', sep="\t") #additional info file
    return freq_df, alleles_df

export_objects()


#barplot
freq,alleles=export_objects()
freq= freq.sort_index()
freq.plot.bar(rot=0,color={"#E6E6FA","#778899"})
plt.savefig('homoGenes_barplot.png', dpi=100)
alleles= alleles.sort_index()
alleles.plot.bar(rot=90,stacked=True,color={"#E6E6FA","#778899"})
plt.tight_layout()
plt.savefig('homoAlleles_barplot.png',pad_inches=4, dpi=300)
plt.cla()

#boxplot
df_a = merge("freq_hla.tsv",re,"resp_a")
df_st= merge("freq_hlaST.tsv",re,"resp_st")
bp_a = df_a.boxplot(column=['diff'])
plt.savefig("diff_a_boxplot.png")
plt.cla()
bp_st = df_st.boxplot(column=['diff'])
plt.savefig("diff_st_boxplot.png")


import pandas as pd
import re
import os
import sys

#----------------------------retrieving data-------------------------------
#amminoacid sequences of alleles
f1=open("sequence.fa","r")
sequences=f1.read()

#grantham matrix
f2="distance.tsv"
matrix= pd.read_csv(f2, sep='\t', index_col=0)

#results from optitype
f3="hla.tsv"
hla= pd.read_csv(f3, sep='\t', header = 0, index_col=0)

#consider couple of alleles for the three genes
def get_alleles(patient):
    L_alleles=[]
    for gene in ['A','B','C']:
        gene1 = re.sub(r'\W+', '', patient[f'{gene}1'])
        gene2 = re.sub(r'\W+', '', patient[f'{gene}2'])
        L_alleles.append((gene1,gene2))
    return L_alleles

#retrieve the sequence for each allele called in get_distance()
def get_seq(allele):
    start=sequences.find(allele)
    end=sequences.find('>',start)
    substring=sequences[start+5:end]
    seq=substring.replace("\n", "")
    return seq

#----------------------------sequence distance---------------------------

#sum the distance between single couples of amminoacids (same position in sequence) and normalize the total distance for the sequence lenght
def get_distance(gene_alleles):
    D=0
    seq1=get_seq(gene_alleles[0])
    seq2=get_seq(gene_alleles[1])
    if len(seq1)!=0 and len(seq2)!=0 and len(seq1)==len(seq2):
        for index in range(len(seq1)):
            aa1=seq1[index]
            aa2=seq2[index]
            d=matrix[aa1][aa2]/len(seq1)
            D+=d
    else:
        D='nan'
    return D

#------------------------------mean distance------------------------------

#calculate mean distance among the three genes
def calc_mean(A,B,C):
    if A=='nan' or B=='nan' or C=='nan':
        mean='nan'
    else:
        mean=(A+B+C)/3
    return mean

#-----------------------distance matrix (sampleXgene)----------------------

def distance_matrix():
    table = []
    for allele, patient in hla.iterrows():
        p_name = patient.name
        alleles = get_alleles(patient)
        A_dist= get_distance(alleles[0])
        B_dist= get_distance(alleles[1])   
        C_dist= get_distance(alleles[2])
        mean_dist = calc_mean(A_dist,B_dist,C_dist)
        line=[p_name,A_dist,B_dist,C_dist,mean_dist]
        table.append(line)
    df=pd.DataFrame(table,columns=['patient','A','B','C','mean'])
    df=df.set_index('patient')
    df.to_csv('HED.tsv', sep='\t')

distance_matrix()    

#----------------------categorical matrix (sampleXgene)--------------------

def assign_category(matrix):


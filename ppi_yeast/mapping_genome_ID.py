import pandas as pd
import numpy as np
import sys,os


## read files
ppi_filename = 'ppi_interactome_goldbi.txt'
id_filename = 'ppi_interactome_goldbi.tab'    #id_file contains 8 fields

with open(ppi_filename,'r') as f1:
    ppi_file = f1.readlines()
    ppi_file = [x.strip() for x in ppi_file]

df = pd.read_csv(id_filename,delimiter='\t')



for row in ppi_file:
    pr1 = row.split('\t')[0]
    pr2 = row.split('\t')[1]

    df1,df2 = df[df['GenomeID']==pr1],df[df['GenomeID']==pr2]
    if (len(df1)!=0) & (len(df2)!=0):
        uniprot_pr1, uniprot_pr2 = df1['Entry'].values[0], df2['Entry'].values[0]

        with open('uniprot_ppi.txt','a') as f:
            f.write(uniprot_pr1+'\t'+uniprot_pr2+'\n')



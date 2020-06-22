#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
import os
import time
import pandas as pd
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Insert amino acids btw individual MSAs to gain flexibility for modelling in trRosetta. ''')
parser.add_argument('--infile', nargs=1, type= str, default=sys.stdin, help = 'Path to infile (MSA in a3m).')
parser.add_argument('--aa_pattern', nargs=1, type= str, default=sys.stdin, help = 'Amino acid pattern to insert.')
parser.add_argument('--pattern_length', nargs=1, type= int, default=sys.stdin, help = 'Length of amino acid pattern to insert.')
parser.add_argument('--seqlens', nargs=1, type= str, default=sys.stdin, help = 'Path to file with sequence lengths for APC.')

###FUNCTIONS###
def flex_msa(infile, l1, aa_pattern, aa_length):
    '''Read a3m MSA and insert aa pattern for flexibility in modelling'''
    flexed = []#Save extracted msa
    with open(infile, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            line = line.rstrip()
            flexed.append(line[:l1]+aa_pattern*aa_length+line[l1:])

    return flexed
		
###MAIN###
args = parser.parse_args()
infile = args.infile[0]
aa_pattern = args.aa_pattern[0]
pattern_length = args.pattern_length[0]
seqlens = pd.read_csv(args.seqlens[0],sep=' ')

#Get protein names
names = infile.split('/')[-1].split('.')[0].split('_')
#Get protein lengths
l1 = seqlens[seqlens['Protein']==names[0]]['Length'].values[0]
l2 = seqlens[seqlens['Protein']==names[1]]['Length'].values[0]
name_pair = names[0]+'_'+names[1]
#Insert flexibility in msa
flexed = flex_msa(infile, l1, aa_pattern, aa_length)
print(flexed)

#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
import os
from collections import Counter
import time
import pandas as pd
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Calculate Mutual Information in MSA. ''')
parser.add_argument('--infile', nargs=1, type= str, default=sys.stdin, help = 'Path to infile (MSA in a3m).')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory.')
parser.add_argument('--seqlens', nargs=1, type= str, default=sys.stdin, help = 'Path to file with sequence lengths for APC.')

###FUNCTIONS###
def read_a3m(infile):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
             'G': 6,'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,'N': 12,
             'O': 21, 'P': 13,'Q': 14, 'R': 15, 'S': 16, 'T': 17,
             'V': 18, 'W': 19, 'Y': 20,'U': 21, 'Z': 21, 'X': 21, 'J': 21}
    max_gap_fraction=0.9
 
    seq = ''
    parsed = []#Save extracted msa
    with open(infile, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if seq != '':
                    gap_fraction = seq.count('-') / float(len(seq))
                    if gap_fraction <= max_gap_fraction:#Only use the lines with less than 90 % gaps
                        parsed.append([mapping.get(ch, 21) for ch in seq if not ch.islower()])
                seq = ''
                continue
            seq += line.rstrip()
        gap_fraction = seq.count('-') / float(len(seq))
        if gap_fraction <= max_gap_fraction:#Only use the lines with less than 90 % gaps
            parsed.append([mapping.get(ch, 21) for ch in seq if not ch.islower()])
    
    print (len(parsed), len(parsed[0]))

    return np.array(parsed, dtype=np.int8, order='F')

def calc_mi(a3m_matrix, name_pair, outdir):
    '''Calculate the MI of a MSA in a3m format
    '''
    n,m = a3m_matrix.shape #n rows(sequences) m columns(amino acids)
    mi_matrix = np.zeros((m,m))

    #Calculate all individual frequencies
    ind_freq = calc_pi(a3m_matrix,m)
    #Calculate the MI
    #Time process
    t1 = time.clock()
    for i in range(m):#All amino acids
        for j in range(i+1,m):#upper triangular columns
            #MI
            Sij = MI(a3m_matrix,i,j,ind_freq)
            mi_matrix[i,j]=Sij#upper triangular assignment
            mi_matrix[j,i]=Sij#lower triangular assignment

    return mi_matrix

def calc_pi(sequences,m):
    '''Calculate the freq of observing aa x in position i for all i.
    m = number of amino acids'''
    ind_freq = []
    for i in range(m):
        ind_freq.append(Counter(s[i] for s in sequences))
    return ind_freq

def MI(sequences,i,j,ind_freq):
    '''Calculate MI'''
    sequences = [s for s in sequences if not 21 in [s[i],s[j]]]
    Pi = ind_freq[i] #The precalculated frequency of all amino acids in position i
    Pj = ind_freq[j] #The precalculated frequency of all amino acids in position j
    Pij = Counter((s[i],s[j]) for s in sequences)
    #log10 or ln?
    return sum(Pij[(x,y)]*np.log(Pij[(x,y)]/(Pi[x]*Pj[y])) for x,y in Pij)

#MAIN
args = parser.parse_args()
infile = args.infile[0]
outdir = args.outdir[0]

names = infile.split('/')[-1].split('.')[0].split('_')
name_pair = names[0]+'_'+names[1]
#Read in msa
a3m_matrix = read_a3m(infile)
#Calculate mututal information
mi_matrix = calc_mi(a3m_matrix, name_pair, outdir)
np.save('array_py.npy', mi_matrix)


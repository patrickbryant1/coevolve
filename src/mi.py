#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
import os
from collections import Counter
import matplotlib.pyplot as plt
import time
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Calculate Mutual Information in MSA. ''')
parser.add_argument('--infile', nargs=1, type= str, default=sys.stdin, help = 'Path to infile (MSA in a3m).')
parser.add_argument('--outfile', nargs=1, type= str, default=sys.stdin, help = 'Path to output filename.')

###FUNCTIONS###
def read_a3m(infile):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
             'G': 6,'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,'N': 12,
             'O': 21, 'P': 13,'Q': 14, 'R': 15, 'S': 16, 'T': 17,
             'V': 18, 'W': 19, 'Y': 20,'U': 21, 'Z': 21, 'X': 21, 'J': 21}
    max_gap_fraction=0.9
    parsed = []#Save extracted msa
    with open(infile, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            line = line.rstrip()
            gap_fraction = line.count('-') / float(len(line))
            if gap_fraction <= max_gap_fraction:#Only use the lines with less than 90 % gaps
                parsed.append([mapping.get(ch, 22) for ch in line if not ch.islower()])

    return np.array(parsed, dtype=np.int8, order='F')

def calc_mi(a3m_matrix):
    '''Calculate the MI of a MSA in a3m format
    '''
    n,m = a3m_matrix.shape #n rows(sequences) m columns(amino acids)
    mi_matrix = np.zeros((m,m))
    #Time process
    t1 = time.clock()
    for i in range(m):#All amino acids
        for j in range(i,m):#upper triangular columns
            #MI
            Sij = MI(a3m_matrix,i,j)
            mi_matrix[i,j]=Sij
    t2 = time.clock()
    print('MI calculation took ', t2-t1, ' s')
    pdb.set_trace()
    plt.imshow(mi_matrix)

    return None

def MI(sequences,i,j):
    sequences = [s for s in sequences if not '-' in [s[i],s[j]]]
    Pi = Counter(s[i] for s in sequences)
    Pj = Counter(s[j] for s in sequences)
    Pij = Counter((s[i],s[j]) for s in sequences)
    return sum(Pij[(x,y)]*np.log(Pij[(x,y)]/(Pi[x]*Pj[y])) for x,y in Pij)

#MAIN
args = parser.parse_args()
infile = args.infile[0]
outfile = args.outfile[0]
#Read in msa
a3m_matrix = read_a3m(infile)
calc_mi(a3m_matrix)
pdb.set_trace()

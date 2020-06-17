#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
import os
from collections import Counter
import matplotlib.pyplot as plt
import time
import pandas as pd
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Perform protein level APC. ''')
parser.add_argument('--infile', nargs=1, type= str, default=sys.stdin, help = 'Path to infile (MSA in a3m).')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory.')
parser.add_argument('--seqlens', nargs=1, type= str, default=sys.stdin, help = 'Path to file with sequence lengths for APC.')



def protein_level_APC(mi_matrix, l1, l2, name_pair, outdir):
    '''
    l1 = length of protein 1
    l2 = length of protein 2
    For  both  MI-and  global-statistic-based  methods,  certain  positions  tend  to  show
    higher  co-evolutionary  signal  with  all  other  positions  due  to  intrinsic  properties
    such  as  higher  entropy.  Average ProductCorrection (APC) was used to solve such problems
    by subtracting the “average product” from the originalmeasurement Sij of co-evolutionary signal:
    '''

    mi_matrix = np.load(outdir+'O13297_Q01159.npy', allow_pickle=True)
    
    corrected_mi = np.zeros(mi_matrix.shape)
    def sum_and_correct(matrix_sel):
        Si_sum = np.sum(matrix_sel,axis=0)#sum of signal over rows
        Sj_sum = np.sum(matrix_sel,axis=1)#sum of signal over columns
        Sij_sum = np.sum(matrix_sel) #sum over complete area
        apc = -np.outer(Si_sum,Sj_sum)/Sij_sum
        return apc

    #Intra protein corrections (aa within the same protein)
    #Protein 1
    matrix_sel = mi_matrix[:l1,:l1]
    apc = sum_and_correct(matrix_sel)
    corrected_mi[:l1,:l1] = mi_matrix[:l1,:l1]-apc
    #Protein 2
    matrix_sel = mi_matrix[l2:,l2:]
    apc = sum_and_correct(matrix_sel)
    corrected_mi[l2:,l2:] = mi_matrix[l2:,l2:]-apc
    #Inter protein corrections (aa btw proteins)
    #Lower triangular
    matrix_sel = mi_matrix[l1:,:l1]
    apc = sum_and_correct(matrix_sel)
    corrected_mi[l1:,:l1] = mi_matrix[l1:,:l1]-apc.T
    #Upper triangular
    matrix_sel = mi_matrix[:l1,l1:]
    apc = sum_and_correct(matrix_sel)
    corrected_mi[:l1,l1:] = mi_matrix[:l1,l1:]-apc.T
    #Save matrix
    np.save(outdir+name_pair+'_per_protein_corr.npy',corrected_mi)
    #Plot
    plot_mi(corrected_mi, outdir+name_pair+'_per_protein_corr.png')
    #Get average MI for later, further APC
    print('Average MI,'+name_pair+','+str(np.average(corrected_mi)))

#MAIN
args = parser.parse_args()
infile = args.infile[0]
outdir = args.outdir[0]
seqlens = pd.read_csv(args.seqlens[0],sep=' ')
#Get protein names
names = infile.split('/')[-1].split('.')[0].split('_')
#Get protein lengths
l1 = seqlens[seqlens['Protein']==names[0]]['Length'].values[0]
l2 = seqlens[seqlens['Protein']==names[1]]['Length'].values[0]
name_pair = names[0]+'_'+names[1]

#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
import os
import matplotlib.pyplot as plt
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Calculate Mutual Information in MSA. ''')
parser.add_argument('--infile1', nargs=1, type= str, default=sys.stdin, help = 'Path to MI matrix 1 in .npy.')
parser.add_argument('--infile2', nargs=1, type= str, default=sys.stdin, help = 'Path to MI matrix 2 in .npy.')
parser.add_argument('--outname', nargs=1, type= str, default=sys.stdin, help = 'Outname, including path to outdir.')

###FUNCTIONS###
def plot_diff(matrix, outname):
    '''Visualize matrix'''
    fig, ax = plt.subplots(figsize=(12,12))
    plt.imshow(matrix)
    cbar = plt.colorbar()
    fig.savefig(outname, format='png', dpi=300)
    plt.close()



#MAIN
args = parser.parse_args()
infile1 = args.infile1[0]
infile2 = args.infile2[0]
outname= args.outname[0]
#Load matrix1
try:
    matrix1 = np.load(infile1, allow_pickle=True)
except:
    matrix1 = np.loadtxt(infile1)
#Load matrix2
try:
    matrix2 = np.load(infile2, allow_pickle=True)
except:
    matrix2 = np.loadtxt(infile2)
pdb.set_trace()
plot_diff(matrix1-matrix2, outname)

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
parser.add_argument('--infile', nargs=1, type= str, default=sys.stdin, help = 'Path to MI matrix in .npy.')
parser.add_argument('--outname', nargs=1, type= str, default=sys.stdin, help = 'Outname, including path to outdir.')

###FUNCTIONS###
def plot_mi(matrix, outname):
    '''Visualize matrix'''
    fig, ax = plt.subplots(figsize=(12,12))
    plt.imshow(matrix)
    cbar = plt.colorbar()
    fig.savefig(outname, format='png', dpi=300)
    plt.close()



#MAIN
args = parser.parse_args()
infile = args.infile[0]
outname= args.outname[0]

matrix = np.load(infile, allow_pickle=True)
plot_mi(matrix, outname)

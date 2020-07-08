#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Visuaize results from model using google mobility data and most of the ICL response team model''')

parser.add_argument('--interacting_pairs', nargs=1, type= str, default=sys.stdin, help = 'Path to interacting pairs (Y2H_union).')
parser.add_argument('--nf5_pairs', nargs=1, type= str, default=sys.stdin, help = 'Path to pairs in paired alignments with Nf5.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

def read_interacting(interacting_id1, interacting_id2, all_ids_in_interactions, nf5_pairs):
    '''Read all pairs with n5 above or equal to 5
    '''
    true_pairs_in_n5 = []
    fetched_pairs = [] #Fetch the pairs that have one of the ids in the true interacting pairs
    with open(nf5_pairs, 'r') as file:
        for line in file:
            line = line.split()
            id1 = line[0]
            id2 = line[1]
            if id1 in
            pdb.set_trace()

#####MAIN#####
#Set font size
matplotlib.rcParams.update({'font.size': 7})
args = parser.parse_args()
interacting_pairs = pd.read_csv(args.interacting_pairs[0], sep = '\t', header = None)
interacting_id1= np.array(interacting_pairs[0])
interacting_id2= np.array(interacting_pairs[1])
all_ids_in_interactions =  np.concatenate([interacting_id1,interacting_id2])
nf5_pairs = args.nf5_pairs[0]
outdir = args.outdir[0]

read_interacting(interacting_id1, interacting_id2, all_ids_in_interactions, nf5_pairs)

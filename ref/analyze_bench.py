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
parser.add_argument('--id_conversion', nargs=1, type= str, default=sys.stdin, help = 'Tab delimited file with ID conversion.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

def format_id_conversion(id_conversion):
    '''Split all gene names and relate to protein id
    '''
    all_conversion = {}
    for i in range(len(id_conversion)):
        row = id_conversion.loc[i]
        protein = row['Entry']
        gene_names = row['Gene names']
        for name in gene_names.split():
            all_conversion[name]=protein

    return all_conversion


def covert_gene_to_protein_id(interacting_ids,id_conversion):
    '''Convert the gene ids to protein ids
    '''
    converted_ids = []
    for id in interacting_ids:
        try:
            converted_ids.append(id_conversion[id])
        except:
            converted_ids.append('missing')
            print(id)
    print('Fraction successfully converted:',len(converted_ids)/len(interacting_ids))

    return converted_ids

def read_interacting(interacting_id1, interacting_id2, all_ids_in_interactions, nf5_pairs):
    '''Read all pairs with n5 above or equal to 5
    '''
    true_pairs_in_n5 = []
    fetched_pairs = [] #Fetch the pairs that have one of the ids in the true interacting pairs
    for id1 in interacting_id1:
        print(id1)
        if id1 in np.array(nf5_pairs[0]):
            pdb.set_trace()
        if id1 in np.array(nf5_pairs[1]):
            pdb.set_trace()
    pdb.set_trace()
    for id2 in interacting_id2:
        if id2 in np.array(nf5_pairs[0]):
            pdb.set_trace()
        if id2 in np.array(nf5_pairs[1]):
            pdb.set_trace()

    pdb.set_trace()
#####MAIN#####
#Set font size
matplotlib.rcParams.update({'font.size': 7})
args = parser.parse_args()
interacting_pairs = pd.read_csv(args.interacting_pairs[0], sep = '\t', header = None)
interacting_id1= np.array(interacting_pairs[0])
interacting_id2= np.array(interacting_pairs[1])
all_ids_in_interactions =  np.unique(np.concatenate([interacting_id1,interacting_id2]))
nf5_pairs = pd.read_csv(args.nf5_pairs[0], sep = ' ', header = None)
id_conversion = pd.read_csv(args.id_conversion[0],sep = '\t')
outdir = args.outdir[0]

#Get all gene names for the id conversion
id_conversion = format_id_conversion(id_conversion)
#Convert all gene ids to protein ids
converted_ids1 = covert_gene_to_protein_id(interacting_id1,id_conversion)
converted_ids2 = covert_gene_to_protein_id(interacting_id2,id_conversion)
pdb.set_trace()
#Map joined msa ids to positive set
read_interacting(id_conversion, interacting_id1, interacting_id2, all_ids_in_interactions, nf5_pairs)

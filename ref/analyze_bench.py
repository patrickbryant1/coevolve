#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np
import seaborn as sns


import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Visuaize results from model using google mobility data and most of the ICL response team model''')

parser.add_argument('--interacting_pairs', nargs=1, type= str, default=sys.stdin, help = 'Path to interacting pairs.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

def read_and_format_data(datadir, countries, days_to_simulate, end_date, covariate_names):


#####MAIN#####
#Set font size
matplotlib.rcParams.update({'font.size': 7})
args = parser.parse_args()
interacting_pairs = args.interacting_pairs[0]
outdir = args.outdir[0]

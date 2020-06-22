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
parser.add_argument('--aa_length', nargs=1, type= int, default=sys.stdin, help = 'Length of amino acid pattern to insert.')
parser.add_argument('--seqlens', nargs=1, type= str, default=sys.stdin, help = 'Path to file with sequence lengths for APC.')



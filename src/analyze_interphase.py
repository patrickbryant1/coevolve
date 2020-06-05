#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Analyze the signal in docked interphases. ''')
parser.add_argument('--docked_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with docked complexes and their pdb files.')
parser.add_argument('--signal_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with coevolution signals such as MI.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


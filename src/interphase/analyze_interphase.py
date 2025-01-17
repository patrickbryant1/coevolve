#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import argparse
import sys
import os
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Analyze the signal in docked interphases. ''')
parser.add_argument('--docked_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with docked complexes and their pdb files.')
parser.add_argument('--signal_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with coevolution signals such as MI.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


##FUNCTIONS###
def parse_atm_record(line):
	'''Parse line in pdb file'''
	record = defaultdict()
	record['name'] = line[0:6].strip()
	record['atm_no'] = int(line[6:11])
	record['atm_name'] = line[12:16].strip()
	record['res_name'] = line[17:20].strip()
	record['chain'] = line[21]
	record['res_no'] = int(line[22:26])
	record['insert'] = line[26].strip()
	record['resid'] = line[22:29]
	record['x'] = float(line[30:38])
	record['y'] = float(line[38:46])
	record['z'] = float(line[46:54])
	record['occ'] = float(line[54:60])
	record['B'] = float(line[60:66])

	return record

def read_pdb(pdbfile):
	'''Read pdb file
	'''

	three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'UNK': 'X'}
	sequence = ''
	res_numbers = []
	prev_res = -1 #Keep track of potential alternative residues
	with open(pdbfile, 'r') as file:
		for line in file:
			try:
				record = parse_atm_record(line)#Parse line
			except:
				print(line)
				continue
			if record['atm_name'] == 'CA': #If CA
				if record['res_no'] == prev_res: #If prev, skip
					continue
				else: #Save residue otherwise
					prev_res = record['res_no']
					sequence += three_to_one[record['res_name']]
					res_numbers.append(record['res_no'])


	return sequence, res_numbers


def get_interphase_score(docked_dir, signal_dir, outdir):
	'''Calculate the score in the docked interphase
	'''

	docked_ids = ['P32914_P27692_2EXU_AA']

	for di in docked_ids:
		seq1, res1 = read_pdb(docked_dir+di+'_u1.pdb')#protein 1
		int1, int_res1 = read_pdb(docked_dir+di+'_u1_M_i1.pdb') #interphase of protein 1
		pdb.set_trace()



#MAIN
args = parser.parse_args()
docked_dir = args.docked_dir[0]
signal_dir = args.signal_dir[0]
outdir = args.outdir[0]
get_interphase_score(docked_dir, signal_dir, outdir)

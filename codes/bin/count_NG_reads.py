#!/usr/bin/env python3
###
# usage: # python3 count_NG_reads.py DICTIONARY_FILE < input_sequences.stgmixd > results
###

import sys
sys.path.append('../lib')
from barcodelib import inputs
from alignment import seq_distance

if __name__ == '__main__':
	dictionary_file = sys.argv[1]
	stg_dict = dict()
	with open(dictionary_file) as f:
		f.readline()
		for line in inputs(f):
			seq, result = line.split()
			if result == 'OK':
				stg_dict[seq] = True
			else:
				stg_dict[seq] = False
	OK = NG = 0
	input()
	for line in inputs():
		_, r, _, stg = line.split()[:4]
		r = int(r)
		if stg not in stg_dict:
			pass
		elif stg_dict[stg]:
			OK += r
		else:
			NG += r
	total = OK + NG
	print('OK reads:', OK, '({}%)'.format(OK*100/total), sep='\t')
	print('NG reads:', NG, '({}%)'.format(NG*100/total), sep='\t')

#!/usr/bin/env python3
###
# usage: # python3 merge_readerror.py [-n] [-w warnings.output] [-r reference.uniq] < input.uniq > output.uniq
###
import os
import sys
import argparse
from collections import defaultdict
from tqdm import tqdm
import time
import shutil
import itertools
sys.path.append('../lib')
from barcodelib import *

max_errors = 2
expandN_bound = 2
save_interval = 900 # 15min

test_flag = False
if test_flag:
    barcode0 = "TTT"
    barcode1 = "TAT"
    print("distance:", hamming_distance(barcode0, barcode1))
    n = 0
    for neighbor in hamming_neighbors(barcode0,2,min_distance=2):
        n += 1
        print(neighbor)
    print(n)
    exit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--warningout', default=None,
                        help='file to output warnings.')
    parser.add_argument('-n', '--noprogress', action='store_true',
                        help='do not show a progress bar.')
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('--hamming', action='store_true',
                        help='Use Hamming distance, instead of Levenshtein.')
    parser.add_argument('--max_errors', type=int, default=max_errors,
                        help='set max_errors. (default:{})'.format(max_errors))
    args = parser.parse_args()
    if args.warningout:
        warningout = open(args.warningout, 'w')
    else:
        warningout = sys.stderr
    errors = 0
    max_errors = args.max_errors
    if args.hamming:
        get_distance = hamming_distance
        get_neighbors = hamming_neighbors
    else:
        get_distance = levenshtein_distance
        get_neighbors = levenshtein_neighbors

    data = []
    if not args.noheader:
        header = input()
    else:
        header = None
    for line in inputs():
        num, arr = line.split()
        num = int(num)
        if 'N' in arr: # skipN
            continue
        else:
            data.append((num, arr))
    data.sort()
    N = len(data)

    remaining_barcodes = {arr for num, arr in data}
    output = []

    start_i = 0
    if args.noprogress:
        barcode_loop = range(start_i,N)
    else:
        barcode_loop = tqdm(range(start_i,N), total=N, initial=start_i)
    for i in barcode_loop:
        readnum, barcode = data[i]
        remaining_barcodes.remove(barcode)
        for d in get_neighbors(barcode, max_errors):
            if d in remaining_barcodes:
                print('Sequence', barcode, 'is an error of', d, sep='\t', file=warningout)
                break
        else:
            output.append((readnum, barcode))

    if header:
        print(header)
    for readnum, barcode in sorted(output, reverse=True):
        print(readnum, barcode, sep='\t')

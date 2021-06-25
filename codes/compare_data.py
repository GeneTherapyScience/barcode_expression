#!/user/bin/env python3
from barcodelib import *
from random import randrange
from collections import Counter, defaultdict
from tqdm import tqdm, trange
from bisect import bisect_right
import argparse
import scipy
import scipy.stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'file1',
        help='input file1.',
    )
    parser.add_argument(
        'file2',
        help='input file2.',
    )
    parser.add_argument(
        '-r', '--reference',
        help='give reference file.',
    )
    parser.add_argument(
        '-i', '--intersection',
        action='store_true',
        help='take intersection of file1 & file2.',
    )
    args = parser.parse_args()

    files = [args.file1, args.file2]
    reads = []
    for i in range(2):
        with open(files[i]) as f:
            data = inputdata(f)[1]
        reads.append(defaultdict(int))
        for record in data:
            barcode, r = record[:2]
            if r > 0:
                reads[-1][barcode] = r
    if args.intersection:
        barcodes_set = set(reads[0].keys()) & set(reads[1].keys())
    else:
        barcodes_set = set(reads[0].keys()) | set(reads[1].keys())
    reference = defaultdict(int)
    if args.reference:
        with open(args.reference) as f:
            data = inputdata(f)[1]
        for record in data:
            if record[1]>0:
                reference[record[0]] = record[1]
    else:
        for barcode in barcodes_set:
            reference[barcode] = 1

    X = [[] for _ in range(2)]
    for barcode in reference.keys() & barcodes_set:
        for i in range(2):
            X[i].append(reads[i][barcode]/reference[barcode])

    P = scipy.stats.pearsonr(*X)
    S = scipy.stats.spearmanr(*X)
    print('Pearson corr:', P)
    print('Spearman corr:', S)

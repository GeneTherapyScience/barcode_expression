#!/usr/bin/env python3
import sys
sys.path.append('../lib')
from barcodelib import *
from tqdm import tqdm, trange
from collections import defaultdict
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    args = parser.parse_args()

    WSN = 'WSN'*10
    header, data = inputdata(has_header=(not args.noheader))
    N = len(data)
    variations = defaultdict(int)
    readnums = defaultdict(int)
    for i in trange(N):
        barcode, readnum, altered = data[i]
        d = levenshtein_distance(barcode, WSN, None)
        variations[d] += 1
        readnums[d] += readnum

    for dic, label in [(variations, 'Variation'), (readnums, 'Readnum')]:
        dist = sorted([(d,n) for d, n in dic.items()])
        av = sum(d*n for d,n in dist)/sum(n for d, n in dist)
        print('{}-base: avarage = {}, distribution = {}'.format(label, av, dist))

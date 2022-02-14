#!/usr/bin/env python3
import sys
sys.path.append('../lib')
from barcodelib import *
import numpy as np
import argparse
from bisect import bisect_left

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('readnumber', type=int,
                        help='give read numbers.')
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    args = parser.parse_args()

    header, data = inputdata(has_header=(not args.noheader))
    N = len(data)
    S = [0]
    for i in range(N):
        S.append(S[-1]+data[i][1])
    order = np.random.permutation(S[-1])[:args.readnumber] + 1
    sample_reads = [0]*N
    sample_mutations = [0]*N
    for k in order:
        i = bisect_left(S, k) - 1
        m = k - S[i]
        sample_reads[i] += 1
        if m <= data[i][2]:
            sample_mutations[i] += 1

    if header:
        print(header)
    for i in range(N):
        if sample_reads[i] > 0:
            print(
                data[i][0],
                sample_reads[i],
                sample_mutations[i],
                sample_mutations[i]/sample_reads[i],
                sep = '\t',
            )

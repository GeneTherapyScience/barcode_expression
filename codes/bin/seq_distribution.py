#!/user/bin/env python3
from barcodelib import *
from alignment import *
from random import randrange
from collections import Counter, defaultdict
from tqdm import tqdm, trange
import argparse

Nsampling = 10
Npair=0
distance_bound = 30

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nsampling', type=int, default=Nsampling,
                        help='give Nsampling.')
    parser.add_argument('-p', '--pairs', type=int, default=Npair,
                        help='give Nsampling.')
    parser.add_argument('-t', '--template', default=None,
                        help='give template to compare with.')
    parser.add_argument('-d', '--detail', default=None,
                        help='give file to write detail.')
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('--all', action='store_true',
                        help='Take all possible pairs.')
    parser.add_argument('--hamming', action='store_true',
                        help='use Hamming distance, rather than Levenshtein.')
    args = parser.parse_args()
    Nsampling = args.nsampling
    Npair = args.pairs
    if args.hamming:
        distance = hamming_distance
    else:
        distance = levenshtein_distance

    header, data = inputdata(has_header=(not args.noheader))
    N = len(data)
    if args.template:
        count = defaultdict(int)
        detail = dict()
        for i in trange(N):
            seq, readnum = data[i][:2]
            if seq in detail:
                detail[seq][2] += readnum
            else:
                align, mut = seq_alignment(args.template, seq)
                detail[seq] = [align, mut, readnum]
            count[mut] += readnum
        print('score', 'reads', sep='\t')
        for k in sorted(count.keys()):
            print(k, count[k], sep='\t')
        if args.detail:
            with open(args.detail, 'w') as f:
                print('template:', args.template, file=f)
                print('sequence', 'alignment', 'score', 'reads', sep='\t', file=f)
                for k, v in sorted(detail.items(), key=lambda x: (x[1][1],-x[1][2],x[0])):
                    print(k, *v, sep='\t', file=f)
    elif args.all:
        count = defaultdict(int)
        M = N*(N-1)//2
        i = j = 0
        for _ in trange(M):
            j += 1
            if j >= N:
                i += 1
                j = i + 1
            count[distance(data[i][0], data[j][0])] += 1
        print(sorted(count.items()), sep='\t')
    elif Npair > 0:
        count = defaultdict(int)
        for _ in trange(Npair):
            while True:
                m, n = randrange(N), randrange(N)
                if m != n:
                    break
            count[distance(data[m][0], data[n][0])] += 1
        print(sorted(count.items()), sep='\t')
    else:
        results = []
        interest_samples = [randrange(N) for _ in range(Nsampling)]
        for i in trange(Nsampling):
            m = interest_samples[i]
            barcode, readnum, _ = data[m]
            count = defaultdict(int)
            for n in range(N):
                if n != m:
                    count[distance(data[n][0], barcode)] += 1
            results.append((barcode, count))

        for barcode, count in results:
            print(barcode, sorted(count.items()), sep='\t')

#!/usr/bin/env python3
###
# usage: # python3 merge_readerror.py [-n] [-w warnings.output] [-r reference.data] < input.data > output.data
###
import os
import sys
import argparse
from collections import defaultdict
from tqdm import tqdm
import time
import shutil
import itertools
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
    parser.add_argument('-r', '--reference', default=None,
                        help='give reference barcodes.')
    parser.add_argument('--refcolumn', type=int, default=0,
                        help='barcode column in reference (default=0).')
    parser.add_argument('--refgroup', default=None,
                        help='give reference groups.')
    parser.add_argument('-w', '--warningout', default=None,
                        help='file to output warnings.')
    parser.add_argument('-u', '--union', action='store_true',
                        help='pack all relatives by union-find algorithm.')
    parser.add_argument('-g', '--groupout', default=None,
                        help='file to output groups of merged relatives.')
    parser.add_argument('-n', '--noprogress', action='store_true',
                        help='do not show a progress bar.')
    parser.add_argument('-l', '--loadfile', default=None,
                        help='file to load halfway results.')
    parser.add_argument('-m', '--milestonefile', default=None,
                        help='file to save halfway results.')
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('--refnoheader', action='store_true',
                        help='the reference does not include header.')
    parser.add_argument('--skipN', action='store_true',
                        help='skip N-including barcodes.')
    parser.add_argument('--hamming', action='store_true',
                        help='Use Hamming distance, instead of Levenshtein.')
    parser.add_argument('--errors', action='store_true',
                        help='output number of errors.')
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

    header, data = inputdata(has_header=(not args.noheader))
    if args.skipN:
        data = list(filter(lambda record: (not 'N' in record[0]), data))
    N = len(data)

    start_i, halfway_barcodes, merged_readnum, merged_mutations = load_halfway(args.loadfile)
    merged_barcodes = readref(args.reference, column=args.refcolumn, has_header=(not args.refnoheader)) | halfway_barcodes
    uf = UnionFind()
    if args.refgroup:
        with open(args.refgroup) as f:
            for line in inputs(f):
                g = line.split()
                p = g[0]
                merged_barcodes.add(p)
                for c in g[1:]:
                    uf.merge(p, c)
                    merged_barcodes.add(c)

    start_t = int(time.time())
    save_t = start_t + save_interval
    if args.noprogress:
        barcode_loop = range(start_i,N)
    else:
        barcode_loop = tqdm(range(start_i,N), total=N, initial=start_i)
    for i in barcode_loop:
        if time.time() > save_t:
            save_t += save_interval
            save_halfway(args.milestonefile, merged_barcodes, merged_readnum, merged_mutations)

        barcode, readnum, mutations = data[i]
        NN = barcode.count('N')
        if NN <= expandN_bound:
            hit = False
            for c in get_neighbors(barcode, max_errors):
                if NN:
                    loop = N_candidates(c)
                else:
                    loop = [c]
                for d in loop:
                    if d in merged_barcodes:
                        if barcode != d and (not uf.connected(d, barcode)):
                            print('Merge', barcode, 'as', d, sep='\t', file=warningout)
                            uf.merge(d, barcode)
                        if not hit:
                            td = uf.root(d) if args.union else d
                            merged_readnum[td] += readnum
                            merged_mutations[td] += mutations
                            if (not args.reference) and args.union:
                                merged_barcodes.add(barcode)
                            hit = True
                        if args.errors:
                            errors += get_distance(d,barcode,max_errors) * readnum
                        if not args.union:
                            break
                else:
                    continue
                break
            if not hit:
                if args.reference:
                    print('Not found in the reference:', barcode, sep='\t', file=warningout)
                else:
                    if 'N' in barcode:
                        print('N-including barcode', barcode, 'has no parent array.', file=warningout)
                    merged_barcodes.add(barcode)
                    merged_readnum[barcode] += readnum
                    merged_mutations[barcode] += mutations
        else:
            for target in merged_barcodes:
                distance = get_distance(barcode, target, max_errors)
                if distance <= max_errors:
                    merged_readnum[target] += readnum
                    merged_mutations[target] += mutations
                    if args.errors:
                        errors += distance * readnum
                    break
            else:
                if args.reference:
                    print('The sequence {} was not found in the reference.'.format(barcode), file=warningout)
                else:
                    print('N-including barcode', barcode, 'has no parent array.', file=warningout)
                    merged_barcodes.add(barcode)
                    merged_readnum[barcode] += readnum
                    merged_mutations[barcode] += mutations

    if header:
        print(header)
    if args.union:
        for barcode in merged_barcodes:
            root = uf.root(barcode)
            if root != barcode:
                print('Merge', barcode, 'as', root, sep='\t', file=warningout)
                merged_readnum[root] += merged_readnum[barcode]
                merged_mutations[root] += merged_mutations[barcode]
                del merged_readnum[barcode], merged_mutations[barcode]
    for barcode in sorted(merged_readnum.keys()):
        readnum, mutations = merged_readnum[barcode], merged_mutations[barcode]
        if readnum == 0:
            ratio = 'NA'
        elif mutations == 0:
            ratio = 0
        elif mutations == readnum:
            ratio = 1
        else:
            ratio = mutations / readnum
        print(barcode, readnum, mutations, ratio, sep='\t')

    if args.groupout:
        groups = defaultdict(list)
        for barcode in uf.elements():
            root = uf.root(barcode)
            if root == barcode:
                groups[root]
            else:
                groups[root].append(barcode)
        with open(args.groupout, 'w') as gout:
            for p, c in groups.items():
                print(p, *c, sep='\t', file=gout)

    if args.errors:
        print(file=sys.stderr)
        print('Read-errors:', errors, file=sys.stderr)

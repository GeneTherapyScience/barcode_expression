#!/usr/bin/python3
###
# usage: # python3 merge_readerror.py [--reference reference.data] < input.data > output.data
###
import sys
import argparse
from itertools import product
from collections import defaultdict

input = lambda: sys.stdin.readline().rstrip()
ATGC = sorted(['A', 'T', 'G', 'C'])
CGTA = ATGC[::-1]

def levenshtein_1(barcode):
    N = len(barcode)
    done = {barcode}
    yield barcode
    for i in range(N):
        yield barcode[:i] + barcode[i+1:]
        for n in ATGC:
            if n != barcode[i]:
                c = barcode[:i] + n + barcode[i+1:]
                if not c in done:
                    done.add(c)
                    yield c
    for i in range(N+1):
        for n in ATGC:
            c = barcode[:i] + n + barcode[i:]
            if not c in done:
                done.add(c)
                yield c

def levenshtein_upto2(barcode):
    done = set()
    for c in levenshtein_1(barcode):
        done.add(c)
        yield c
    for c in levenshtein_1(barcode):
        for d in levenshtein_2(c):
            if not d in done:
                done.add(d)
                yield d

def N_candidates(barcode):
    k = barcode.count('N')
    if k == 0:
        yield barcode
    else:
        barcode_list = list(barcode)
        Npos = [i for i, x in enumerate(barcode) if x == 'N']
        pool = [(n, 0) for n in CGTA]
        while pool:
            c, i = pool.pop()
            barcode_list[Npos[i]] = c
            if i == k-1:
                yield str(barcode_list)
            else:
                for n in CGTA:
                    pool.append((n, i+1))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', default=None,
                        help='give reference barcodes.')
    args = parser.parse_args()

    header = input()
    data = []
    while True:
        try:
            line = input().strip()
            if line[0] == '#' or len(line)==0:
                continue
            else:
                barcode, readnum, altered = line.split()[:3]
                readnum = 0 if readnum=='NA' else int(readnum)
                altered = 0 if altered=='NA' else int(altered)
                data.append((barcode, readnum, altered))
        except EOFError:
            break
    data.sort(key=lambda x: x[0].count('N'), -x[1], x[0])

    merged_barcodes = set()
    original_barcode = dict()
    merged_readnum = defaultdict(int)
    merged_altered = defaultdict(int)
    if args.reference:
        if os.path.exists(args.reference):
            with open(args.reference) as f:
                while True:
                    try:
                        line = input().strip()
                        if line[0] == '#' or len(line)==0:
                            continue
                        else:
                            barcode = line.split()[0]
                            merged_barcodes.add(barcode)
                            for d in levenshtein_upto2(barcode):
                                original_barcode[d] = barcode
                    except EOFError:
                        break
        else:
            print('Error: reference file "{}" does not exist.'.format(args.reference),
                  file=sys.stderr)
    else:
        for barcode, readnum, altered in data:
            if not 'N' in barcode:
                if barcode in original_barcode:
                    merged_readnum[original_barcode[barcode]] += readnum
                    merged_altered[original_barcode[barcode]] += altered
                elif args.reference:
                    print('The sequence {} was not found in the reference.'.format(barcode), file=sys.stderr)
                else:
                    merged_barcodes.add(barcode)
                    merged_readnum[barcode] += readnum
                    merged_altered[barcode] += altered
                    original_barcode[barcode] = barcode
                    for d in levenshtein_upto2(barcode):
                        original_barcode[d] = barcode
            else:
                for c in N_candidates(barcode):
                    if c in original_barcode:
                        merged_readnum[original_barcode[c]] += readnum
                        merged_altered[original_barcode[c]] += altered
                        break
                else:
                    print('N-including barcode', barcode, 'has no parent array.', file=sys.stderr)

    for barcode in sorted(merged_barcodes):
        readnum, altered = merged_readnum[barcode], merged_altered[barcode]
        if readnum == 0:
            ratio = 'NA'
        elif altered == 0:
            ratio = 0
        elif altered == readnum:
            ratio = 1
        else:
            ratio = altered / readnum
        print(barcode, readnum, altered, ratio, sep='\t')

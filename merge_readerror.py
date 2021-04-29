#!/usr/bin/python3
###
# usage: # python3 merge_readerror.py [-n] [-w warnings.output] [-r reference.data] < input.data > output.data
###
import os
import sys
import argparse
from collections import defaultdict
from tqdm import trange

max_errors = 2
expandN_bound = 2

ATGC = sorted(['A', 'T', 'G', 'C'])
CGTA = ATGC[::-1]

def levenshtein_neighbors(barcode, distance, inserts=ATGC):
    if distance == 0:
        yield barcode
    else:
        prev = list()
        for b in levenshtein_neighbors(barcode, distance-1):
            prev.append(b)
            yield b
        done = set(prev)
        for b in prev:
            N = len(b)
            for i in range(N):
                yield b[:i] + b[i+1:]
                if b[i] == 'N':
                    continue
                for n in inserts:
                    if n != b[i]:
                        c = b[:i] + n + b[i+1:]
                        if not c in done:
                            done.add(c)
                            yield c
            for i in range(N+1):
                for n in inserts:
                    c = b[:i] + n + b[i:]
                    if not c in done:
                        done.add(c)
                        yield c

def levenshtein_distance(barcode0, barcode1, bound=None):
    M, N = len(barcode0), len(barcode1)
    if bound is None:
        bound = M + N
    if M < N:
        M, N, barcode0, barcode1 = N, M, barcode1, barcode0
    cur = list(range(1,N+1)) + [0]
    for m in range(M):
        k = min(cur)
        if k > bound:
            return k
        prev = cur
        cur = [0]*(N+1)
        cur[-1] = m+1
        print(prev, cur)
        for n in range(N):
            cur[n] = min(
                prev[n-1] + int(barcode0[m]!='N' and barcode1[n]!='N' and barcode0[m]!=barcode1[n]),
                cur[n-1] + 1,
                prev[n] + 1,
                )
    return cur[N-1]

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

def readref(referencefile):
    ret = set()
    with open(referencefile) as f:
        f.readline() # drop header
        while True:
            line = f.readline()
            if not line: # EOF
                break
            line = line.strip()
            if line[0] == '#' or len(line)==0:
                continue
            else:
                barcode = line.split()[0]
                ret.add(barcode)
    return ret

def readstdin():
    while True:
        try:
            line = input().strip()
            if len(line)==0 or line[0] == '#':
                continue
            else:
                yield line
        except EOFError:
            break

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', default=None,
                        help='give reference barcodes.')
    parser.add_argument('-w', '--warningout', default=None,
                        help='file to output warnings.')
    parser.add_argument('-n', '--noprogress', action='store_true',
                        help='do not show a progress bar.')
    args = parser.parse_args()
    if args.warningout:
        warningout = open(args.warningout, 'w')
    else:
        warningout = sys.stderr
    if args.noprogress:
        mrange = range
    else:
        mrange = trange

    header = input()
    data = []
    for line in readstdin():
        barcode, readnum, altered = line.split()[:3]
        readnum = 0 if readnum=='NA' else int(readnum)
        altered = 0 if altered=='NA' else int(altered)
        data.append((barcode, readnum, altered))
    data.sort(key=lambda x: (x[0].count('N'), -x[1], x[0]))
    N = len(data)

    if args.reference:
        try:
            merged_barcodes = readref(args.reference)
        except Exception as e:
            print(e, file=warningout)
            exit()
    else:
        merged_barcodes = set()

    merged_readnum = defaultdict(int)
    merged_altered = defaultdict(int)
    for i in mrange(N):
        barcode, readnum, altered = data[i]
        if barcode.count('N') <= expandN_bound:
            for c in levenshtein_neighbors(barcode, max_errors, inserts=['N']):
                for d in N_candidates(c):
                    if d in merged_barcodes:
                        merged_readnum[d] += readnum
                        merged_altered[d] += altered
                        break
                else:
                    continue
                break
            else:
                if args.reference:
                    print('The sequence {} was not found in the reference.'.format(barcode), file=warningout)
                else:
                    if 'N' in barcode:
                        print('N-including barcode', barcode, 'has no parent array.', file=warningout)
                    merged_barcodes.add(barcode)
                    merged_readnum[barcode] += readnum
                    merged_altered[barcode] += altered
        else:
            for target in merged_barcodes:
                if levenshtein_distance(barcode, target, max_errors) <= max_errors:
                    merged_readnum[target] += readnum
                    merged_altered[target] += altered
                    break
            else:
                if args.reference:
                    print('The sequence {} was not found in the reference.'.format(barcode), file=warningout)
                else:
                    print('N-including barcode', barcode, 'has no parent array.', file=warningout)
                    merged_barcodes.add(barcode)
                    merged_readnum[barcode] += readnum
                    merged_altered[barcode] += altered

    print(header)
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

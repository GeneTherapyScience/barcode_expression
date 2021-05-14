#!/usr/bin/python3
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

max_errors = 2
expandN_bound = 2
save_interval = 900 # 15min

ATGC = sorted('ATGC')
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
                c = b[:i] + b[i+1:]
                if not c in done:
                    done.add(c)
                    yield c
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

DNA_match_pairs = set()
for b in 'ATGC':
    DNA_match_pairs.add(b+b)
for n in 'ATGCWSN':
    DNA_match_pairs |= {'N'+n, n+'N'}
for w in 'ATW':
    DNA_match_pairs |= {'W'+w, w+'W'}
for s in 'GCS':
    DNA_match_pairs |= {'S'+s, s+'S'}

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
        for n in range(N):
            cur[n] = min(
                prev[n-1] + int(barcode0[m]+barcode1[n] not in DNA_match_pairs),
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

def inputs(f=sys.stdin):
    while True:
        line = f.readline()
        if not line: # EOF
            break
        else:
            yield line

def readref(referencefile):
    ret = set()
    if referencefile:
        with open(referencefile) as f:
            f.readline() # drop header
            for line in inputs(f):
                line = line.strip()
                if line[0] == '#' or len(line)==0:
                    continue
                else:
                    barcode = line.split()[0]
                    ret.add(barcode)
    return ret

def inputdata(f=sys.stdin, has_header=True):
    if has_header:
        header = f.readline().rstrip()
    else:
        header = None
    data = []
    for line in inputs(f):
        barcode, readnum, mutations = line.split()[:3]
        readnum = 0 if readnum=='NA' else int(readnum)
        mutations = 0 if mutations=='NA' else int(mutations)
        data.append((barcode, readnum, mutations))
    data.sort(key=lambda x: (x[0].count('N'), -x[1], x[0]))
    return header, data

def load_halfway(loadfile):
    halfway_barcodes = set()
    halfway_readnum = defaultdict(int)
    halfway_mutations = defaultdict(int)
    start_i = 0
    if loadfile:
        with open(loadfile) as f:
            start_i = int(f.readline())
            for line in inputs(f):
                barcode, readnum, mutations = line.split()
                halfway_barcodes.add(barcode)
                halfway_readnum[barcode] = int(readnum)
                halfway_mutations[barcode] = int(mutations)
    return start_i, halfway_barcodes, halfway_readnum, halfway_mutations

def save_halfway(savefile, merged_barcodes, merged_readnum, merged_mutations):
    if not savefile:
        return False
    tmpfile = savefile+'.tmp'
    with open(tmpfile, 'w') as f:
        print(i, file=f)
        for barcode in sorted(merged_barcodes):
            readnum, mutations = merged_readnum[barcode], merged_mutations[barcode]
            print(barcode, readnum, mutations, sep='\t', file=f)
    shutil.move(tmpfile, savefile)
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', default=None,
                        help='give reference barcodes.')
    parser.add_argument('-w', '--warningout', default=None,
                        help='file to output warnings.')
    parser.add_argument('-n', '--noprogress', action='store_true',
                        help='do not show a progress bar.')
    parser.add_argument('-l', '--loadfile', default=None,
                        help='file to load halfway results.')
    parser.add_argument('-m', '--milestonefile', default=None,
                        help='file to save halfway results.')
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('--skipN', action='store_true',
                        help='skip N-including barcodes.')
    args = parser.parse_args()
    if args.warningout:
        warningout = open(args.warningout, 'a')
    else:
        warningout = sys.stderr

    header, data = inputdata(has_header=(not args.noheader))
    if args.skipN:
        data = list(filter(lambda record: (not 'N' in record[0]), data))
    N = len(data)

    start_i, halfway_barcodes, merged_readnum, merged_mutations = load_halfway(args.loadfile)
    merged_barcodes = readref(args.reference) | halfway_barcodes

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
            for c in levenshtein_neighbors(barcode, max_errors):
                if NN:
                    loop = N_candidates(c)
                else:
                    loop = [c]
                for d in loop:
                    if d in merged_barcodes:
                        merged_readnum[d] += readnum
                        merged_mutations[d] += mutations
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
                    merged_mutations[barcode] += mutations
        else:
            for target in merged_barcodes:
                if levenshtein_distance(barcode, target, max_errors) <= max_errors:
                    merged_readnum[target] += readnum
                    merged_mutations[target] += mutations
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
    for barcode in sorted(merged_barcodes):
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

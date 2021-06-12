#!/user/bin/env python3
from merge_readerror import *
from random import randrange
from collections import Counter
from tqdm import tqdm, trange
from bisect import bisect_right
import argparse

spieces = 10**5
cells = 2*(10**6)

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--reference', default=None)
parser.add_argument('-s', '--spieces', type=int, default=spieces)
parser.add_argument('-c', '--cells', type=int, default=cells)
parser.add_argument('-i', '--incubate', type=int, default=0)
args = parser.parse_args()
reference = args.reference
spieces = args.spieces
cells = args.cells

if reference:
    with open(reference) as f:
        header, data = inputdata(f)
        spieces = len(data)
        S = [record[1] for record in data]
        M = sum(S)
else:
    S = [0]*spieces
    for i in range(cells):
        S[randrange(spieces)] += 1
    M = cells
Acc = [0]*(spieces)
for i in range(spieces):
    Acc[i] = Acc[i-1] + S[i]

for t in range(args.incubate):
    T = [0]*spieces
    for i in range(cells):
        r = randrange(M)
        p = bisect_right(Acc, r)
        T[p] += 1
    S = T
    Acc = [0]*spieces
    for s in range(spieces):
        Acc[s] = Acc[s-1] + S[s]
    M = cells
cells = sum(S)

R = [0]*spieces
var = 0
reads = [10**5, 3*(10**5), 10**6, 3*(10**6), 65*(10**5), 
10**7, 2*(10**7),
3*(10**7), 10**8
]
reads_set = set(reads)
if args.incubate:
    spieces = len(list(filter(lambda x: x>0, S)))
    print(spieces, "variations after incubation.")
for i in range(1,reads[-1]+1):
    if i in reads_set:
        print('{} reads\t{} variations'.format(i, var))
    r = randrange(cells)
    p = bisect_right(Acc, r)
    if R[p] == 0:
        var += 1
        if var == spieces:
            print('{} reads\t{} variations'.format(i,var))
            break
    R[p] += 1

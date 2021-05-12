#!/user/bin/python3
from merge_readerror import *
from random import randrange
from collections import Counter
from tqdm import tqdm, trange
from bisect import bisect_right
import argparse

spieces = 10**5
cells = 2*(10**6)

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--spieces', type=int, default=spieces)
parser.add_argument('-c', '--cells', type=int, default=cells)
args = parser.parse_args()
spieces = args.spieces
cells = args.cells

S = [0]*spieces
for i in range(cells):
    S[randrange(spieces)] += 1
Acc = [0]*(spieces)
for i in range(spieces):
    Acc[i] = Acc[i-1] + S[i]

R = [0]*spieces
var = 0
reads = [10**5, 3*(10**5), 10**6, 3*(10**6), 10**7, 3*(10**8), 10**8]
reads_set = set(reads)
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
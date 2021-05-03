#!/usr/bin/python
from merge_readerror import *
from tqdm import tqdm, trange
from collections import defaultdict
# import argparse

WSN = 'WSN'*10
header, data = inputdata()
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

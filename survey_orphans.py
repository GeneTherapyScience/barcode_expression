#!/user/bin/python3
from merge_readerror import *
from random import randrange
from collections import Counter
from tqdm import tqdm, trange

Nsampling = 10
distance_bound = 10

header, data = inputdata()
parents = []
orphans = []
for barcode, readnum, altered in data:
    if readnum == 1 or 'N' in barcode:
        orphans.append((barcode, readnum, altered))
    else:
        parents.append((barcode, readnum, altered))

#interest = orphans
interest = parents
N = len(interest)
interest_samples = [randrange(N) for _ in range(Nsampling)]
results = []
for i in trange(Nsampling):
    barcode, readnum, altered = interest[interest_samples[i]]
    min_d = distance_bound+1
    parent = 'N/A'
    for target, target_readnum, target_altered in parents:
        if barcode == target:
            continue
        d = levenshtein_distance(barcode, target, min_d)
        if d < min_d:
            min_d = d
            parent = target
            # break
    # print(min_d, barcode, parent, sep='\t')
    results.append((min_d, barcode, parent))

for d, barcode, target in results:
    print(d, barcode, target, sep='\t')

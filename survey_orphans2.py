#!/user/bin/python3
from merge_readerror import *
from random import randrange
from collections import Counter
from tqdm import tqdm, trange
import argparse

Nsampling = 10
default_distance_bound = 2

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--distance', type=int, default=default_distance_bound)
parser.add_argument('-n', '--number', type=int, default=Nsampling)
args = parser.parse_args()
Nsampling = args.number
distance_bound = args.distance

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
    n = r = 0
    for target, target_readnum, target_altered in parents:
        d = levenshtein_distance(barcode, target, distance_bound+1)
        if d <= distance_bound:
            n += 1
            r += target_readnum
            # break
    # print(min_d, barcode, parent, sep='\t')
    results.append((barcode, n, r))

print('barcode', '#variations', '#reads', sep='\t')
for barcode, n, r in results:
    print(barcode, n, r, sep='\t')
print('Avarage', *[sum(x[i] for x in results)/Nsampling for i in [1,2]], sep='\t')

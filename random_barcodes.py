#!/usr/bin/python3
import argparse
from random import randrange

barcode_length = 30
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--number', type=int, default=10**5)
parser.add_argument('-w', '--wsn', action='store_true')
args = parser.parse_args()

print('WSN')
def wsn_random(i):
    if i % 3 == 0:
        return 'AT'[randrange(2)]
    elif i % 3 == 1:
        return 'GC'[randrange(2)]
    else:
        return 'ATGC'[randrange(4)]

if args.wsn:
    for _ in range(args.number):
        print(*[ wsn_random(i) for i in range(barcode_length)], '\t2\t0', sep='')
else:
    for _ in range(args.number):
        print(*[ 'ATGC'[randrange(4)] for i in range(barcode_length)], '\t2\t0', sep='')

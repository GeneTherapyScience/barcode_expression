#!/user/bin/env python
import sys
sys.path.append('../lib')
from barcodelib import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('reference',
                        help='reference set.')
    parser.add_argument('inputfile', default=None,
                        help='input file.')
    # parser.add_argument('--uniq', action='store_true',
    #                     help='uniq file input.')
    args = parser.parse_args()

    ref_barcodes = set()
    with open(args.reference) as f:
        if not args.noheader:
            f.readline()
        for line in inputs(f):
            ref_barcodes.add(line.split()[1])
    if args.inputfile:
        infile = open(args.inputfile)
    else:
        infile = sys.stdin
    if not args.noheader:
        infile.readline()
    for line in inputs(infile):
        barcode = line.split()[1]
        if barcode in ref_barcodes:
            print(line)

import argparse
import sys
sys.path.append('../lib')
from barcodelib import inputs
from collections import defaultdict
Nwsn = 10

def is_wsn(barcode, n):
    L = 3*n
    if len(barcode) != L:
        return False
    for i in range(L):
        if i % 3 == 0:
            if barcode[i] not in {'A', 'T', 'W', 'N'}:
                return False
        elif i % 3 == 1:
            if barcode[i] not in {'G', 'C', 'S', 'N'}:
                return False
    else:
        return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('--onlybarcode', action='store_true',
                        help='output only barcodes.')
    parser.add_argument('--column', type=int, default=0,
                        help='column where to read barcode.')
    parser.add_argument('--filelist', action='store_true',
                        help='read filelist from stdin & output overwrap.')
    args = parser.parse_args()

    if args.filelist:
        Nf = 0
        appearance = defaultdict(int)
        for file in inputs():
            Nf += 1
            with open(file) as f:
                if not args.noheader:
                    header = f.readline().strip()
                for line in inputs(f):
                    record = line.split()
                    barcode = record[args.column]
                    if is_wsn(barcode, Nwsn):
                        appearance[barcode] += 1
        for b in sorted([b for b, n in appearance.items() if n == Nf]):
            print(b)
    else:
        if not args.noheader:
            header = input()
            if not args.onlybarcode:
                print(header)
        for line in inputs():
            record = line.split()
            barcode = record[args.column]
            if is_wsn(barcode, Nwsn):
                if args.onlybarcode:
                    print(barcode)
                else:
                    print(*record, sep='\t')

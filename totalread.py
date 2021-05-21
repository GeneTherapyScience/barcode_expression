#!/user/bin/python3
from merge_readerror import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    args = parser.parse_args()

    header, data = inputdata(has_header=(not args.noheader))
    totalread = totalmutation = 0
    totalvariation = 0
    for record in data:
        barcode, read, mutation = record[:3]
        totalread += read
        totalmutation += mutation
        if read > 0:
            totalvariation += 1
    print('variations', 'reads', 'mutation', sep='\t')
    print(totalvariation, totalread, totalmutation, sep='\t')

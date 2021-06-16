#!/user/bin/env python
from merge_readerror import *
from wsnstg2mutinfo import wsn_thres
from bisect import bisect_left
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('-r', '--ratiothres', type=float, default=0,
                        help='give ratio for threshold.')
    args = parser.parse_args()

    header, data = inputdata(has_header=(not args.noheader))
    totalread = totalmutation = 0
    totalvariation = 0
    if args.ratiothres > 0:
        ratiothres_flag = True
        read_arr = []
    else:
        ratiothres_flag = False
    for record in data:
        barcode, read, mutation = record[:3]
        totalread += read
        totalmutation += mutation
        if read > 0:
            totalvariation += 1
        if ratiothres_flag:
            read_arr.append(read)
    titles = ['variations', 'reads', 'mutation']
    values = [totalvariation, totalread, totalmutation]
    if ratiothres_flag:
        read_arr.sort()
        th = wsn_thres(read_arr, ratio=args.ratiothres)
        n = bisect_left(read_arr, th)
        titles.append('limited_variations')
        values.append(totalvariation - n)
    print(*titles, sep='\t')
    print(*values, sep='\t')

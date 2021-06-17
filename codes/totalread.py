#!/user/bin/env python
from merge_readerror import *
from wsnstg2mutinfo import wsn_thres
from bisect import bisect_left
import math
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    parser.add_argument('-r', '--ratiothres', type=float, default=0,
                        help='give ratio for threshold.')
    parser.add_argument('-s', '--shannon', action='store_true',
                        help='show Shannon-variation.')
    args = parser.parse_args()

    header, data = inputdata(has_header=(not args.noheader))
    totalread = totalmutation = 0
    totalvariation = 0
    read_arr = []
    if args.ratiothres > 0:
        ratiothres_flag = True
    else:
        ratiothres_flag = False
    for record in data:
        barcode, read, mutation = record[:3]
        totalread += read
        totalmutation += mutation
        if read > 0:
            totalvariation += 1
        read_arr.append(read)
    titles = ['variations', 'reads', 'mutation']
    values = [totalvariation, totalread, totalmutation]
    if ratiothres_flag:
        read_arr.sort()
        th = wsn_thres(read_arr, ratio=args.ratiothres)
        n = bisect_left(read_arr, th)
        titles.append('limited_variations')
        values.append(totalvariation - n)
    if args.shannon:
        S = -sum(
            ((r/totalread) * math.log(r/totalread) if r>0 else 0)
            for r in read_arr)
        titles.append('Shannon-variation')
        v = math.exp(S)
        if v >= 100:
            v = round(v)
        else:
            v = round(v,1)
        values.append(v)
    print(*titles, sep='\t')
    print(*values, sep='\t')

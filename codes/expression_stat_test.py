from collections import defaultdict
from barcodelib import inputs
import numpy as np
from bisect import bisect_left, bisect_right
from scipy import stats
from math import log10

key = 'insdel_mean'
p_th = 0.05

def get_order(arr, reverse=False):
    order = []
    N = len(arr)
    vals = sorted(arr)
    for v in arr:
        l = bisect_left(vals, v)
        r = bisect_right(vals, v)
        m = (l+r-1)//2 if (l+r)%2 else (l+r-1)/2
        order.append(N-1-m if reverse else m)
    return order

def order2score(x):
    return x

if __name__ == '__main__':
    headers = input().split()
    key_columns = [i for i in range(len(headers)) if headers[i].split(':')[-1]==key]
    W = len(key_columns)
    key_headers = [headers[i] for i in key_columns]
    barcodes = []
    values = [list() for _ in range(W)]
    for line in inputs():
        line_values = line.split()
        barcodes.append(line_values[0])
        for i in range(W):
            values[i].append(line_values[key_columns[i]])
    orders = [get_order(values[i]) for i in range(W)]
    N = len(orders[0])
    p_th_eff = p_th/N
    envs = ['NT', 'DTX', 'Sp']
    K = len(envs)
    env_columns = [list() for _ in range(K)]
    for i in range(W):
        for k in range(K):
            if envs[k] in key_headers[i]:
                env_columns[k].append(i)
                break
    for j in range(N):
        scores = [list() for _ in range(K)]
        for k in range(K):
            for i in env_columns[k]:
                scores[k].append(order2score(orders[i][j]))
        for k in range(K):
            target = scores[k]
            remnant = []
            for l in range(K):
                if l != k:
                    remnant += scores[l]
            t, p = stats.ttest_ind(target, remnant, equal_var = False)
            if p < p_th_eff:
                target_mean = sum(target)/len(target)
                remnant_mean = sum(remnant)/len(remnant)
                if target_mean > remnant_mean:
                    sign = '+'
                else:
                    sign = '-'
                print(barcodes[j], j, k, sign, target_mean, remnant_mean, p*N, sep='\t')

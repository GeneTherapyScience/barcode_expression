from collections import defaultdict
import sys
sys.path.append('../lib')
from barcodelib import inputs
import numpy as np
from bisect import bisect_left, bisect_right
from scipy import stats
import math
import sys

key = 'insdel_reads'

if __name__ == '__main__':
    for filename in sys.argv[1:]:
        with open(filename) as f:
            headers = f.readline().rstrip().split()
            key_columns = [i for i in range(len(headers)) if headers[i].split(':')[-1]==key]
            W = len(key_columns)
            key_headers = [headers[i] for i in key_columns]
            barcodes = []
            values = [list() for _ in range(W)]
            for line in inputs(f):
                line_values = line.split()
                barcodes.append(line_values[0])
                for i in range(W):
                    values[i].append(float(line_values[key_columns[i]]))

            N = len(values[0])
            envs = ['NT', 'DTX', 'Sp']
            K = len(envs)
            env_columns = [list() for _ in range(K)]
            for i in range(W):
                for k in range(K):
                    if envs[k] in key_headers[i]:
                        env_columns[k].append(i)
                        break

            ret = [[0]*(len(env_columns[k])+1) for k in range(K)]
            for j in range(N):
                for k in range(K):
                    n = 0
                    for i in env_columns[k]:
                        if values[i][j]:
                            n += 1
                    ret[k][n] += 1

            for k in range(K):
                print(filename, envs[k], *ret[k], sep='\t')

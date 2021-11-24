from collections import defaultdict
from barcodelib import inputs
import numpy as np
from bisect import bisect_left, bisect_right
from scipy import stats
import math
import sys

key = 'insdel_mean'
p_th = 0.05
STAT_N = 10**8

def get_order(arr, reverse=True):
    order = []
    N = len(arr)
    vals = sorted(arr)
    for v in arr:
        l = bisect_left(vals, v)
        r = bisect_right(vals, v)
        m = (l+r+1)//2 if (l+r)%2 else (l+r+1)/2
        order.append(N+1-m if reverse else m)
    return order

def get_scores(arr):
    scores = []
    N = len(arr)
    vals = sorted(arr)
    M = bisect_right(vals, 0)
    for v in arr:
        if v <= 0:
            s = 0
        else:
            l = bisect_left(vals, v)
            r = bisect_right(vals, v)
            m = (l+r-1)/2
            s = 2*(m+1-M)/(N-M)*(N/(N-M+1))
        scores.append(s)
    return scores

def get_deviation(arr):
    arr = np.array(arr)
    print(arr)
    m = arr.mean()
    sigma = arr.std(ddof=1)
    return (arr - m) / sigma

def order2score(x):
    return x

def test_score(target, remnant):
    return np.log(target).mean() - np.log(remnant).mean()

def reference_distribution(N, Ntarget, Nremnant):
    rng = np.random.Generator(np.random.PCG64())
    ret = [np.log(rng.uniform(0,1,Ntarget)).mean() - np.log(rng.uniform(0,1,Nremnant)).mean() for _ in range(N)]
    ret.sort()
    return ret

class StatJudge:
    def __init__(self, zero_ratio, N=10**8, K=3):
        self.N = N
        alpha = np.array(zero_ratio)
        p = 2/(1-alpha)**2
        rng = np.random.Generator(np.random.PCG64())
        W = len(alpha)
        V = W//3
        self.experience = [list() for _ in range(K)]
        for _ in range(self.N):
            x = (rng.uniform(0,1,W) - alpha) * p
            s0 = np.where(x>0, x, 0)
            s1 = np.array([s0[V*k:V*(k+1)].sum() for k in range(K)])
            S = s1.sum()
            for k in range(K):
                self.experience[k].append(s0[k]*2 - S)
        for k in range(K):
            self.experience[k].sort()

    def __call__(self, s, k):
        return (self.N - bisect_left(self.experience[k], s))/self.N


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
            values[i].append(float(line_values[key_columns[i]]))
    all_scores = [get_scores(values[i]) for i in range(W)]
    zero_ratio = [values[i].count(0)/len(values[i]) for i in range(W)]
    N = len(all_scores[0])
    p_th_eff = p_th/N
    envs = ['NT', 'DTX', 'Sp']
    K = len(envs)
    env_columns = [list() for _ in range(K)]
    for i in range(W):
        for k in range(K):
            if envs[k] in key_headers[i]:
                env_columns[k].append(i)
                break
    # C = len(env_columns[0])
    # reference = reference_distribution(STAT_N-1, C, C*2)
    stat_judge = StatJudge(zero_ratio, N=10**7)
    for j in range(N):
        scores = [list() for _ in range(K)]
        for k in range(K):
            for i in env_columns[k]:
                scores[k].append(all_scores[i][j])
        for k in range(K):
            target = scores[k]
            remnant = []
            for l in range(K):
                if l != k:
                    remnant += scores[l]
            # t, p = stats.ttest_ind(target, remnant, equal_var = False)
            # s = test_score(target, remnant)
            # p = (bisect_left(reference, s) + 1) / STAT_N

#            if p < p_th_eff or 1-p < p_th_eff:
            # target_mean = np.array(target).mean()
            # remnant_mean = np.array(remnant).mean()
            target_mean = np.array(target).mean()
            remnant_mean = np.array(remnant).mean()
            s = target_mean - remnant_mean
            p = stat_judge(s, k)
            print(barcodes[j], j, k, '{:.8f}'.format(target_mean), target, '{:.8f}'.format(remnant_mean), remnant, '{:.8f}'.format(s), p, sep='\t')

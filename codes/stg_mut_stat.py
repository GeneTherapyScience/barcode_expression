import analyze_data01 as ad
from itertools import product
# from scipy.stats import spearmanr, weightedtau
import numpy as np
import pandas as pd
from collections import defaultdict
# from matthijsz.weightedcorr.WeightedCorr import WeightedCorr
from scipy import special
import sys

clones = ['N112', 'P112', 'S98']
env = ['NT', 'DTX', 'Sp']
days = ['D0', 'D7', 'D14']
E, T, K = 3, 3, 3
target_t = 2
t_range = [target_t]
min_reads = 4

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='stat analyses of mut-level in different envs, for each barcode')
    args = parser.parse_args()

    # print('Weighted correlations : Pearson & Spearman')
    for clone in clones:
        print(clone)
        barcodes = set()
        data = ad.get_celllines(clone, datadir='../wsnstg_white40_hamm2B/', common = '.merge.extracted.reformat.white40.filtered.wsnstg')
        for e, t, k in product(range(E), t_range, range(K)):
            barcodes |= set(data[e][t][k].keys())
        result = []
        for barcode in barcodes:
            mu, V, reads = [np.empty((3,3)) for _ in range(3)]
            for e, t, k in product(range(E), t_range, range(K)):
                cur = data[e][t][k][barcode]
                mu[e][k] = cur['insdel_mean']
                V[e][k] = cur['insdel_V']
                reads[e][k] = cur['reads']
            for max_e in range(E):
                # other_e = set(range(E)) - {max_e}
                if min(reads[max_e].sum(),reads.sum()-reads[max_e].sum()) < min_reads:
                    continue
                m0 = (mu[max_e] @ reads[max_e]) / reads[max_e].sum()
                m1 = (np.tensordot(mu,reads,((0,1),(0,1))) - mu[max_e] @ reads[max_e]) / (reads.sum() - reads[max_e].sum())
                V0 = (V[max_e] + mu[max_e]**2) @ reads[max_e] / reads[max_e].sum() - m0*m0
                V1 = (np.tensordot(V+mu**2, reads, ((0,1),(0,1))) - (V[max_e]+mu[max_e]**2) @ reads[max_e]) / (reads.sum() - reads[max_e].sum()) - m1*m1
                # if min(V0,V1) <= 0:
                #     print(V0, V1, file=sys.stderr)
                sigma = np.sqrt(V0 + V1)
                if m0 <= m1:
                    continue
                w = (m0 - m1)/sigma
                result.append((1-special.erf(w), max_e, barcode, reads.sum()))
        result.sort()
        N = 32
        n = min(N, len(result))
        print('Most significant {} are:'.format(n))
        print('p-value', 'env', 'barcode', 'reads', sep='\t')
        for i in range(n):
            print(result[i][0], env[result[i][1]], result[i][2], result[i][3], sep='\t', flush=True)

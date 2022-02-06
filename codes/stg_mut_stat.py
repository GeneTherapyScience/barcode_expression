import analyze_data01 as ad
from itertools import product
# from scipy.stats import spearmanr, weightedtau
import numpy as np
import pandas as pd
from collections import defaultdict
# from matthijsz.weightedcorr.WeightedCorr import WeightedCorr
from scipy import special

clones = ['N112', 'P112', 'S98']
env = ['NT', 'DTX', 'Sp']
days = ['D0', 'D7', 'D14']
E, T, K = 3, 3, 3
target_t = 2
t_range = [target_t]
min_reads = 64

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='stat analyses of mut-level in different envs, for each barcode')
    args = parser.parse_args()

    # print('Weighted correlations : Pearson & Spearman')
    for clone in clones:
        print(clone)
        barcodes = set()
        reads = defaultdict(int)
        data = ad.get_celllines(clone, datadir='../wsnstg_white40_hamm2B/', common = '.merge.extracted.reformat.white40.filtered.wsnstg')
        for e, t, k in product(range(E), t_range, range(K)):
            barcodes |= set(data[e][t][k].keys())
        result = []
        for barcode in barcodes:
            mu, V = np.empty((3,3)), np.empty((3,3))
            for e, t, k in product(range(E), t_range, range(K)):
                cur = data[e][t][k][barcode]
                mu[e][k] = cur['insdel_mean']
                V[e][k] = cur['insdel_V']
                if cur['reads'] < min_reads:
                    continue
            max_e = np.argmax([mu[e][0] for e in range(E)])
            other_e = set(range(E)) - {max_e}
            if min(mu[max_e]) > max([mu[e][k] for e in other_e for k in range(K)]):
                w = min((mu[max_e][k0]-mu[oe][k1])/np.sqrt(V[max_e][k0]+V[oe][k1])
                        for oe, k0, k1 in product(other_e, range(K), range(K)))
                result.append((1-special.erf(w), max_e, barcode))
        result.sort()
        N = 32
        print('{} barcodes, {} are 3-6 structure & all more than {} reads.'.format(len(barcodes),len(result),min_reads))
        n = min(N, len(result))
        print('Most significant {} are:'.format(n))
        print('p-value', 'env', 'barcode', sep='\t')
        for i in range(n):
            print(result[i][0], env[result[i][1]], result[i][2], sep='\t', flush=True)

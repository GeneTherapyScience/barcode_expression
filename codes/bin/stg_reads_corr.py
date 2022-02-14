import sys
sys.path.append('../lib')
sys.path.append('../thirdparty')
import analyze_data01 as ad
from itertools import product
from scipy.stats import spearmanr, weightedtau
import pandas as pd
from collections import defaultdict
from matthijsz.weightedcorr.WeightedCorr import WeightedCorr

clones = ['N112', 'N128', 'P112', 'P114', 'S98', 'S101']
E, T, K = 3, 3, 3
env = ['NT', 'DTX', 'Sp']
days = ['D0', 'D7', 'D14']

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='output reads-insdel correlations')
    parser.add_argument('--table', action='store_true',
                        help='output table')
    parser.add_argument('--pearson', action='store_true',
                        help='pearson correlation')
    parser.add_argument('--spearman', action='store_true',
                        help='spearman correlation')
    args = parser.parse_args()

    # print('Weighted correlations : Pearson & Spearman')
    result = defaultdict(dict)
    for clone in clones:
        data = ad.get_celllines(clone, datadir='../wsnstg_white40_hamm2B/', common = '.merge.extracted.reformat.white40.filtered.wsnstg')
        for e, t, k in product(range(E), range(T), range(K)):
            cur = data[e][t][k]
            reads, insdel_mean = [], []
            for wsn in cur:
                reads.append(cur[wsn]['reads'])
                insdel_mean.append(cur[wsn]['insdel_mean'])
            #corr, pvalue = spearmanr(reads, insdel_mean)
            #corr, pvalue = weightedtau(reads, insdel_mean, weigher=lambda i: weights[i])
            wc = WeightedCorr(x=pd.Series(reads),y=pd.Series(insdel_mean),w=pd.Series(reads))
            corr = []
            if args.pearson:
                corr.append(wc(method='pearson'))
            if args.spearman:
                corr.append(wc(method='spearman'))
            result[clone][(e,t,k)] = corr

    if args.table:
        for clone in clones:
            output = []
            header = [clone] + env + ['total']
            output.append(header)
            c = sum([result[clone][(0,0,k)][0] for k in range(K)])/K
            record = [days[0], c] + ['-']*(E-1) + [c]
            output.append(record)
            for t in range(1,T):
                record = [days[t]]
                for e in range(E):
                    c = sum([result[clone][(e,t,k)][0] for k in range(K)])/K
                    record.append(c)
                record.append(sum(record[1:])/E)
                output.append(record)
            record = ['+'.join(days[1:])]
            for e in range(1,E+2):
                record.append(sum(output[t][e] for t in range(2,T+1))/(T-1))
            output.append(record)
            for record in output:
                print(*record, sep='\t')
            print(flush=True)
    else:
        for e, t, k in product(range(E), range(T), range(K)):
            print(clone, days[t], env[e], k+1, *corr, sep='\t', flush=True)

import sys
sys.path.append('../lib')
import numpy as np
import analyze_data01 as ad
import sys
from collections import defaultdict
from itertools import product

def get_scores(data,e,t,wsns,keys):
    E, T, K = len(data), len(data[0]), len(data[0][0])
    X = np.array([ad.get_subdata_readmut(data,e,t,k,*keys,only_alive=False,wsns=wsns) for k in range(K)])
    return X

def avg_scores(data,e,t,wsns,keys):
    E, T, K = len(data), len(data[0]), len(data[0][0])
    X = get_scores(data,e,t,wsns,['reads']+keys).sum(axis=2)
    for k in range(K):
        X[k,1:] /= X[k,0]
    return X

def compare_stg(data,e,t,wsns,keys=['muts','total_mut']):
    P0 = avg_scores(data,0,t,wsns,keys).mean(axis=0)[1:]
    print(*P0, sep='\t')
    if e != 0:
        P1 = avg_scores(data,e,t,wsns,keys).mean(axis=0)[1:]
        print(*P1, sep='\t')
    else:
        P1 = P0
    return tuple([P1[i]/P0[i] for i in range(len(P0))])

def base_wsns(data):
    day0 = data[0][0]
    wsns = set()
    for data_k in day0:
        wsns |= set(data_k.keys())
    return list(wsns)

def output_scores(
    outfile,
    samples=[('NANOG',112),('NANOG',128),('NANOG',151),('SOX2',98),('SOX2',101),('POU5F1',112),('POU5F1',114)],
    datadir='../wsnstg_white40_hamm2B/',
    ratio=0,
    common='.merge.extracted.reformat.white40.filtered.wsnstg',
    ):
    keys = ['mut', 'ins', 'del', 'mis', 'mut_reads', 'ins_reads', 'del_reads', 'mis_reads']
    env = ['NT','DTX','Sp']
    days = ['Day0','Day7','Day14']
    f = open(outfile, 'w')
    print('sample', 'reads', *keys, sep='\t', file=f)
    for gene, num in samples:
        long_name = '{}_{}'.format(gene,num)
        short_name = '{}{}'.format(gene[0],num)
        tot_data = ad.get_celllines(long_name, datadir=datadir, ratio=ratio, common=common)
        data = ad.get_celllines(short_name, datadir=datadir, ratio=ratio, common=common)
        wsns = base_wsns(data)
        scores = avg_scores(tot_data,0,0,keys=keys,wsns=wsns)
        print(long_name, int(scores[0][0]), *scores[0][1:], sep='\t', file=f)
        scores = avg_scores(data,0,0,keys=keys,wsns=wsns)
        for k in range(3):
            print('{}{}_{}'.format(short_name,days[0],k+1), int(scores[k][0]), *scores[k][1:], sep='\t', file=f)
        data = ad.get_celllines(short_name, datadir=datadir, ratio=ratio, common=common)
        wsns = base_wsns(data)
        for e in range(3):
            for t in range(1,3):
                scores = avg_scores(data,e,t,keys=keys,wsns=wsns)
                for k in range(3):
                    print('{}{}{}_{}'.format(short_name,days[t],env[e],k+1), int(scores[k][0]), *scores[k][1:], sep='\t', file=f)
    f.close()

def output_histograms(
    outfile="stg_{}_histograms.tsv",
    samples=[('NANOG',112),('NANOG',128),('NANOG',151),('SOX2',98),('SOX2',101),('POU5F1',112),('POU5F1',114)],
    datadir='../wsnstg_white40_hamm2B/',
    ratio=0,
    common='.merge.extracted.reformat.white40.filtered.wsnstg',
    ):
    env = ['NT','DTX','Sp']
    days = ['Day0','Day7','Day14']
    divs = ['mut', 'ins', 'del', 'mis']
    M = 80
    f = dict()
    for q in divs:
        f[q] = open(outfile.format(q), 'w')
        print('sample',*[i/2 for i in range(M)], sep='\t', file=f[q])
    for gene, num in samples:
        long_name = '{}_{}'.format(gene,num)
        short_name = '{}{}'.format(gene[0],num)
        tot_hist = ad.get_stg_histogram(long_name, datadir=datadir, ratio=ratio, common=common)
        hist = ad.get_stg_histogram(short_name, datadir=datadir, ratio=ratio, common=common)
        for q in divs:
            h = [0]*M
            for p,n in tot_hist[0][0][0][q].items():
                h[int(2*p)] = n/tot_hist[0][0][0]['reads']
            print(long_name, *h, sep='\t', file=f[q])
            for k in range(3):
                h = [0]*M
                for p,n in hist[0][0][k][q].items():
                    h[int(2*p)] = n/hist[0][0][k]['reads']
                print('{}{}_{}'.format(short_name,days[0],k+1), *h, sep='\t', file=f[q])
            for e in range(3):
                for t in range(1,3):
                    for k in range(3):
                        h = [0]*M
                        for p,n in hist[e][t][k][q].items():
                            h[int(2*p)] = n/hist[e][t][k]['reads']
                        print('{}{}{}_{}'.format(short_name,days[t],env[e],k+1), *h, sep='\t', file=f[q])
    for q in divs:
        f[q].close()

def get_stg_cdf(
    outfile="stg_{}_cdfs.tsv",
    samples=[('NANOG',112),('NANOG',128),('NANOG',151),('SOX2',98),('SOX2',101),('POU5F1',112),('POU5F1',114)],
    datadir='../wsnstg_white40_hamm2B/',
    ratio=0,
    common='.merge.extracted.reformat.white40.filtered.wsnstg',
    ):
    divs = ['mut_mean', 'ins_mean', 'del_mean', 'mis_mean', 'mut_ratio', 'ins_ratio', 'del_ratio', 'mis_ratio']
    env = ['NT','DTX','Sp']
    days = ['Day0','Day7','Day14']
    f = dict()
    for s,q in enumerate(divs):
        f[q] = open(outfile.format(divs[s]), 'w')
        print('sample',*[i/100 for i in range(101)], sep='\t', file=f[q])
    for gene, num in samples:
        long_name = '{}_{}'.format(gene,num)
        short_name = '{}{}'.format(gene[0],num)
        tot_data = ad.get_celllines(long_name, datadir=datadir, ratio=ratio, common=common)
        data = ad.get_celllines(short_name, datadir=datadir, ratio=ratio, common=common)
        wsns = base_wsns(data)
        W = len(wsns)
        scores = get_scores(tot_data,0,0,keys=['reads']+divs,wsns=wsns)
        scores.sort(axis=2)
        for s,q in enumerate(divs):
            print(long_name, *[scores[0][s+1][round((W-1)*i/100)] for i in range(101)], sep='\t', file=f[q])
        scores = get_scores(data,0,0,keys=['reads']+divs,wsns=wsns)
        scores.sort(axis=2)
        for k in range(3):
            for s,q in enumerate(divs):
                print('{}{}_{}'.format(short_name,days[0],k+1), *[scores[k][s+1][round((W-1)*i/100)] for i in range(101)], sep='\t', file=f[q])
        for e in range(3):
            for t in range(1,3):
                scores = get_scores(data,e,t,keys=['reads']+divs,wsns=wsns)
                scores.sort(axis=2)
                for k in range(3):
                    print('{}{}_{}'.format(short_name,days[0],k+1), *[scores[k][s+1][round((W-1)*i/100)] for i in range(101)], sep='\t', file=f[q])
    for q in divs:
        f[q].close()

def make_wsn_mut_table(
    outfile="../wsn_mut_table/wsn_mut_table_{}_a.tsv",
    clones=[('NANOG',112),('NANOG',128),('NANOG',151),('SOX2',98),('SOX2',101),('POU5F1',112),('POU5F1',114)],
    datadir='../wsnstg_white40_hamm2B/',
    ratio=0,
    common='.merge.extracted.reformat.white40.filtered.wsnstg',
    divs = ['reads','insdel_reads','insdel_ratio','insdel','insdel_mean'],
    ):
    # divs = ['mut_mean', 'ins_mean', 'del_mean', 'mis_mean', 'mut_ratio', 'ins_ratio', 'del_ratio', 'mis_ratio']
    E, T, K = 3, 3, 3
    env = ['NT','DTX','Sp']
    days = ['Day0','Day7','Day14']
    for gene, num in clones:
        long_name = '{}_{}'.format(gene,num)
        short_name = '{}{}'.format(gene[0],num)

        header = ([long_name]
        + ['total:{}'.format(div) for div in divs]
        + [
            '{}{}_{}:{}'.format(days[t],env[e] if t > 0 else '',k+1,div)
            for e in range(E) for t in range(T) for k in range(K)
            for div in divs if not (e>0 and t==0)
        ])
        table = defaultdict(list)
        tot_data = ad.get_celllines(long_name, datadir=datadir, ratio=ratio, common=common)
        data = ad.get_celllines(short_name, datadir=datadir, ratio=ratio, common=common)
        wsns = set(tot_data[0][0][0].keys())
        for e, t, k in product(range(E),range(T),range(K)):
            if e>0 and t==0:
                continue
            wsns |= data[e][t][k].keys()
        W = len(wsns)

        for wsn in wsns:
            table[wsn] = (
                [tot_data[0][0][0][wsn][div] for div in divs]
                + [data[e][t][k][wsn][div]
                    for e in range(E) for t in range(T) for k in range(K) 
                    for div in divs if not (e>0 and t==0)])

        with open(outfile.format(long_name), 'w') as f:
            print(*header, sep='\t', file=f)
            for wsn, scores in sorted(table.items(), key=lambda x: x[1], reverse=True):
                print(wsn, *scores, sep='\t', file=f)


if __name__ == '__main__':
    # get_stg_cdf()
    make_wsn_mut_table()

from scipy.special import comb
import itertools
import numpy as np
from collections import defaultdict

def start_end_reads(
    data,
    e, t,
):
    data0 = data[0][0] # Day0
    reads0 = defaultdict(int)
    for k in range(3):
        for wsn, d in data0[k].items():
            reads0[wsn] += d['reads']
    wsns0 = sorted(reads0.keys(), key=lambda x: (reads0[x], x))
    start_reads = [reads0[wsn] for wsn in wsns0]
    end_reads = [[data[e][t][k][wsn]['reads'] for wsn in wsns0] for k in range(3)]
    return start_reads, end_reads

def wsn_updown(
    start_reads,
    end_reads,
    Ndiv = 8,
    threshold = 0, # percentage
    display = True,
):
    M = len(start_reads)
    rng = np.random.Generator(np.random.PCG64()) # 同数の場合ランダムに
    start_order = (np.array(start_reads)+rng.uniform(0,1/2,M)).argsort()
    end_order = (np.array(end_reads)+rng.uniform(0,1/2,(3,M))).argsort(axis=1)
    start_pos = [-1]*M
    end_pos = [[-1]*M for _ in range(3)]
    for p in range(M):
        w = start_order[p]
        start_pos[w] = p
        for k in range(3):
            w = end_order[k][p]
            end_pos[k][w] = p
    up = [0]*M
    down = [0]*M
    for w in range(M):
        sp = start_pos[w]
        ep = [end_pos[k][w] for k in range(3)]
        if (np.array(ep)-sp > M*threshold).sum()==3:
            up[w] = 1
        elif (np.array(ep)-sp < -M*threshold).sum()==3:
            down[w] = 1
    U, D = sum(up), sum(down)
    if display:
        print("total:")
        print(U, D, U/M, D/M, (U+D)/M, sep='\t')
        print("{}-classes:".format(Ndiv))
    d = np.linspace(0,M,Ndiv+1)
    ret = []
    for i in range(len(d)-1):
        prev, cur = map(int,d[i:i+2])
        U, D = sum(up[prev:cur]), sum(down[prev:cur])
        K = cur-prev
        ret.append((U, D, U/K, D/K, (U+D)/K))
        if display:
            print(*ret[-1], sep='\t')
    return ret

def start_end_reads2(
    data,
    e, t,
):
    data0 = data[0][0] # Day0
    reads0 = defaultdict(int)
    for k in range(3):
        for wsn, d in data0[k].items():
            reads0[wsn] += d['reads']
    wsns0 = sorted(reads0.keys(), key=lambda x: (reads0[x], x))
    start_reads = [[data0[k][wsn]['reads'] for wsn in wsns0] for k in range(3)]
    end_reads = [[data[e][t][k][wsn]['reads'] for wsn in wsns0] for k in range(3)]
    return start_reads, end_reads

def wsn_updown2(
    start_reads,
    end_reads,
    Ndiv = 8,
    threshold = 0, # percentage
    display = True,
):
    M = len(start_reads[0])
    start_k = len(start_reads)
    end_k = len(end_reads)
    rng = np.random.Generator(np.random.PCG64()) # 同数の場合ランダムに
    start_order = (np.array(start_reads)+rng.uniform(0,1/2,(start_k,M))).argsort(axis=1)
    end_order = (np.array(end_reads)+rng.uniform(0,1/2,(end_k,M))).argsort(axis=1)
    start_pos = [[-1]*M for _ in range(start_k)]
    end_pos = [[-1]*M for _ in range(end_k)]
    for p in range(M):
        for k in range(start_k):
            w = start_order[k][p]
            start_pos[k][w] = p
        for k in range(end_k):
            w = end_order[k][p]
            end_pos[k][w] = p
    up = [0]*M
    down = [0]*M
    for w in range(M):
        sp_max = max(start_pos[k][w] for k in range(start_k))
        sp_min = min(start_pos[k][w] for k in range(start_k))
        ep_max = max(end_pos[k][w] for k in range(end_k))
        ep_min = min(end_pos[k][w] for k in range(end_k))

        if ep_min-sp_max > M*threshold:
            up[w] = 1
        elif ep_max-sp_min < -M*threshold:
            down[w] = 1
    d = np.linspace(0,M,Ndiv+1)
    ret = []
    if display:
        print("{}-classes:".format(Ndiv))

    for i in range(len(d)-1):
        prev, cur = map(int,d[i:i+2])
        U, D = sum(up[prev:cur]), sum(down[prev:cur])
        K = cur-prev
        ret.append((K, U, D, U/K, D/K, (U+D)/K))
        if display:
            print(*ret[-1], sep='\t')

    U, D = sum(up), sum(down)
    ret.append((M, U, D, U/M, D/M, (U+D)/M))
    if display:
        print("total:")
        print(M, U, D, U/M, D/M, (U+D)/M, sep='\t')

    return ret

def wsn_updown3(
    start_reads,
    end_reads,
    display=True,
):
    M = len(start_reads[0])
    start_k = len(start_reads)
    end_k = len(end_reads)
    rng = np.random.Generator(np.random.PCG64()) # 同数の場合ランダムに
    start_order = (np.array(start_reads)+rng.uniform(0,1/2,(start_k,M))).argsort(axis=1)
    end_order = (np.array(end_reads)+rng.uniform(0,1/2,(end_k,M))).argsort(axis=1)
    start_pos = -np.ones((start_k,M))
    end_pos = -np.ones((end_k,M))
    for p in range(M):
        for k in range(start_k):
            w = start_order[k][p]
            start_pos[k][w] = p
        for k in range(end_k):
            w = end_order[k][p]
            end_pos[k][w] = p
    start_pos = start_pos + rng.uniform(0,1/2,(start_k, M))
    end_pos = end_pos + rng.uniform(0,1/2,(end_k, M))
    N = start_k+end_k
    L = comb(N, start_k, exact=True)
    dist = [0]*L
    c2n = defaultdict(int)
    for i, c in enumerate(itertools.combinations(range(N),start_k)):
        c2n[c] = i
    score = 0
    for w in range(M):
        order = np.array([start_pos[k][w] for k in range(start_k)]
                          +[end_pos[k][w] for k in range(end_k)]
                         ).argsort()
        cmb = tuple(np.where(order < start_k)[0])
        dist[c2n[cmb]] += 1
        score += abs(15 - 2*sum(cmb))
    control = [M/L]*L
    score /= M

    if display:
        print([m/M for m in dist])
        print("score =", score)
    
    ret = scipy.stats.chisquare(dist, control)

    return ret


import sys
import argparse
import re
import analyze_data01 as ad
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--datadir', default='./wsnstg_white40_hamm2B/',
                        help='set datadir.')
    parser.add_argument('--common', default='.merge.extracted.reformat.white40.filtered.wsnstg',
                        help='common suffix for filenames.')
    parser.add_argument('--dictfile', default='./stginfo/whitelist.sorted.stgmixd',
                        help='stgmixd dict-file.')
    # parser.add_argument('--perread', action='store_true',
    #                     help='output counts per read.')
    args = parser.parse_args()

    clones = [
        "N112", "N128", "N151", "P112", "P114", "S98", "S101",
    ]
    env = ['NT', 'DTX', 'Sp']
    time = ['D0', 'D7', 'D14']

    print('clone', 'day', 'env', 'total-var', 'all-up', 'all-down', 'all-up-ratio', 'all-down-ratio', 'all-(up/down)-ratio', sep='\t')
    for clone in clones:
        data = ad.get_celllines(clone, datadir=args.datadir, dictfile=args.dictfile, ratio=10**(-5), common = args.common)
        for e in range(3):
            for t in range(1,3):
                ret=wsn_updown2(*start_end_reads2(data,e,t),display=False,Ndiv=1)[-1]
                print(clone, time[t], env[e], *ret, sep='\t')

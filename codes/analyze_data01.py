import os
from collections import defaultdict
from barcodelib import inputs
import numpy as np

def get_filelist(sample, datadir='../wsnstg_white40/', common = '.merge.extracted.reformat.white40.wsnstg'):
    days = ['Day0', 'D7', 'D14']
    env = ['NT', 'DTX', 'Sp']
    filelist = []
    for e in range(3):
        filelist.append([])
        for t in range(3):
            filelist[-1].append([])
            sample_str = sample
            if e == 1 and t > 0: # DTX
                if sample_str[1] == '_':
                    sample_str = sample_str[0] + sample_str[2:]
            else:
                if sample_str[1] != '_':
                    sample_str = sample_str[0] + '_' + sample_str[1:]

            for k in range(3):
                if t == 0:
                    filename = os.path.join(datadir, '{}_{}_{}{}'.format(days[t], sample_str, k+1, common))
                else:
                    filename = os.path.join(datadir, '{}_{}_{}_{}{}'.format(days[t], sample_str, env[e], k+1, common))
                filelist[-1][-1].append(filename)
    return filelist

from itertools import product
from wsnstg2mutinfo import wsn_thres
from barcodelib import get_mixd_dictionary, mixd_distance

def get_celllines(sample, datadir='../wsnstg_white40/', dictfile='../stginfo/whitelist.sorted.stgmixd', ratio=10**(-5), common = '.merge.extracted.reformat.white40.wsnstg'):
    files = get_filelist(sample, datadir, common)
    mixd_dictionary = get_mixd_dictionary(dictfile)
    result = defaultdict(lambda: [[],[],[]])
    for e, t, k in product(range(3), repeat=3):
        filename = files[e][t][k]
        wsn_data = defaultdict(lambda: defaultdict(int))
        wsn_dist = defaultdict(lambda: defaultdict(int))
        with open(filename) as f:
            f.readline()
            for line in inputs(f):
                wsn, stg, r = line.split()
                r = int(r)
                mixd = mixd_dictionary[stg]
                d, di, dx, dd = mixd_distance(mixd)
                wsn_data[wsn]['reads'] += r
                wsn_data[wsn]['muts'] += r * int(d > 0)
                wsn_data[wsn]['ins'] += di*r
                wsn_data[wsn]['del'] += dd*r
                wsn_data[wsn]['mis'] += dx*r
                wsn_dist[wsn][d] += r
        th = wsn_thres([wsn_data[wsn]['reads'] for wsn in wsn_data.keys()], ratio)
        for wsn in wsn_dist.keys():
            wsn_data[wsn]['dist'] = wsn_dist[wsn]
            wsn_data[wsn]['total_mut'] = sum(d*r for d,r in wsn_data[wsn]['dist'].items())
            wsn_data[wsn]['mean_all'] = wsn_data[wsn]['total_mut']/wsn_data[wsn]['reads']
            wsn_data[wsn]['mean_mut'] = wsn_data[wsn]['total_mut']/wsn_data[wsn]['muts'] if wsn_data[wsn]['muts']>0 else None
            wsn_data[wsn]['ins_mean'] = wsn_data[wsn]['ins']/wsn_data[wsn]['reads']
            wsn_data[wsn]['del_mean'] = wsn_data[wsn]['del']/wsn_data[wsn]['reads']
            wsn_data[wsn]['mis_mean'] = wsn_data[wsn]['mis']/wsn_data[wsn]['reads']
            wsn_data[wsn]['alive'] = True if wsn_data[wsn]['reads'] >= th else False

        result[e][t].append(wsn_data)
    return result

def get_subdata_readmut(data, e, t, k, *keys, only_alive=True, wsns=None):
    subdata = data[e][t]
    if isinstance(k,int):
        K = [k]
    else:
        K = k
    subdata_readmut = [
        [subdata[k][wsn][key] for key in keys]
        for k in K for wsn in (wsns or subdata[k].keys())
        if (((not only_alive) or subdata[k][wsn]['alive'])
            and (None not in [subdata[k][wsn][key] for key in keys]))
    ]
    X = [[] for _ in range(len(keys))]
    for xy in subdata_readmut:
        for i in range(len(keys)):
            X[i].append(xy[i])
    X = np.array(X)
    return X

if __name__ == '__main__':
    print(get_filelist('N112'))

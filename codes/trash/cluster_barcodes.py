from barcodelib import inputs
import numpy as np
from sklearn.cluster import KMeans
import sys
import argparse

if __name__ == '__main__':
    seed = 0
    K = 5
    key = 'insdel_mean'

    print('#', sys.argv)

    parser = argparse.ArgumentParser()
    parser.add_argument('-K', '--clusters', type=int, default=K, help='specify cluster number (default: {}).'.format(K))
    parser.add_argument('--key', default=key)
    parser.add_argument('inputfile', type=str, help='input tsv file')
    args = parser.parse_args()
    K = args.clusters
    key = args.key
    inputfile = args.inputfile

    f = open(inputfile)
    headers = f.readline().split()
    key_columns = [i for i in range(len(headers)) if headers[i].split(':')[-1]==key][1:]
    reads_columns = [i for i in range(len(headers)) if headers[i].split(':')[-1]=='reads'][1:]
    samplenames = [headers[c].split(':')[0] for c in key_columns]
    W = len(key_columns)
    key_headers = [headers[i] for i in key_columns]
    barcodes = []
    values = []
    reads = []
    for line in inputs(f):
        line_values = line.split()
        barcodes.append(line_values[0])
        values.append([float(line_values[c]) for c in key_columns])
        reads.append(sum(int(line_values[c]) for c in reads_columns))
    # X = np.array(values)
    # weights = np.array(reads)

    km = KMeans(n_clusters=K, random_state=seed)
    km.fit(values, sample_weight=reads)
    print(*samplenames, '\t')
    for y in (np.log10(km.cluster_centers_ + 1e-10)*100).astype(int)/100:
        print(*y, sep='\t')

    # envs = ['NT', 'DTX', 'Sp']
    # K = len(envs)
    # env_columns = [list() for _ in range(K)]
    # for i in range(W):
    #     for k in range(K):
    #         if envs[k] in key_headers[i]:
    #             env_columns[k].append(i)
    #             break


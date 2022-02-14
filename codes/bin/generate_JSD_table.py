import sys
sys.path.append('../lib')
from scipy.spatial.distance import jensenshannon
from barcodelib import inputs
import numpy as np

top_num = 10
envs = ['NT', 'DTX', 'Sp']
E, K = 3, 3

if __name__ == '__main__':
    reference_wsnstg = sys.argv[1] # common wsn-file
    datafile = sys.argv[2]
    output_total = sys.argv[3]
    output_top = sys.argv[4]
    valid_barcodes = set()
    with open(reference_wsnstg) as f:
        for line in inputs(f):
            barcode = line.split()[0]
            if barcode[0] == '#':
                continue
            else:
                valid_barcodes.add(barcode)

    with open(datafile) as f:
        while True:
            record = f.readline().split()
            if record[0] != '#':
                header = record
                break

        clone_name = header[0]
        sample_names = (
            ['Day0_{}'.format(k+1) for k in range(K)] + 
            ['Day{}{}_{}'.format(n,e,k+1)
                for n in [7,14] for e in envs for k in range(K)]
        )
        column_names_set = {'{}:reads'.format(s) for s in sample_names }
        reads_columns = [c for c in range(len(header)) if header[c] in column_names_set]
        M = len(reads_columns)
        reads = [list() for _ in range(M)]

        for line in inputs(f):
            record = line.split()
            barcode = record[0]
            if barcode in valid_barcodes:
                for m in range(M):
                    reads[m].append(int(record[reads_columns[m]]))

    top_barcodes = set(np.argsort(reads)[:,-top_num:].flatten())
    top_reads = np.array(reads)[:,np.array(list(top_barcodes))]

    output_files = [output_total, output_top]
    output_data = [reads, top_reads]

    for u in range(2):
        with open(output_files[u], 'w') as f:
            print('#', *sys.argv,
                sep='\t', file=f)
            print(clone_name, *sample_names,
                sep='\t', file=f)
            for h in range(M):
                print(sample_names[h], *[jensenshannon(output_data[u][h],output_data[u][w]) for w in range(M)],
                    sep='\t', file=f)

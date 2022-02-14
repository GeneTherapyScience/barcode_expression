import sys
sys.path.append('../lib')
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

        barcodes = []
        for line in inputs(f):
            record = line.split()
            barcode = record[0]
            if barcode in valid_barcodes:
                for m in range(M):
                    reads[m].append(int(record[reads_columns[m]]))
                barcodes.append(barcode)

    reads = np.array(reads)
    top_barcodes = set(np.argsort(reads)[:,-top_num:].flatten())
    top_order = np.array(sorted(list(top_barcodes), key=lambda i: -reads[:,i].sum()))
    top_reads = np.array(reads)[:,top_order]
    top_ratios = top_reads / top_reads.sum(axis=1)[:,np.newaxis]

    total_order = np.array(sorted(range(len(barcodes)), key=lambda i: -reads[:,i].sum()))
    total_reads = np.array(reads[:,total_order])
    total_ratios = total_reads / total_reads.sum(axis=1)[:,np.newaxis]

    output_files = [output_total, output_top]
    output_data = [
        (total_order, total_reads, total_ratios),
        (top_order, top_reads, top_ratios),
    ]

    column_names = [sample+sp for sample in sample_names for sp in [':reads',':ratio']]
    for u in range(2):
        with open(output_files[u], 'w') as f:
            print('#', *sys.argv,
                sep='\t', file=f)
            print(clone_name, *column_names,
                sep='\t', file=f)
            order, reads, ratios = output_data[u]
            for i in range(len(order)):
                print(barcodes[order[i]], *[reads[w,i] if s==0 else ratios[w,i] for w in range(M) for s in range(2)],
                    sep='\t', file=f)

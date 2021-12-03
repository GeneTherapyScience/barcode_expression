import sys
from barcodelib import inputs
from collections import defaultdict

if __name__ == '__main__':
    print('#', *sys.argv)
    reference_wsnstg = sys.argv[1] # Plasmid wsnstg-file
    datafiles = sys.argv[2:]
    B2N = defaultdict(int)
    with open(reference_wsnstg) as f:
        f.readline()
        for line in inputs(f):
            barcode = line.split()[0]
            B2N[barcode] = 1

    for filename in datafiles:
        with open(filename) as f:
            header = f.readline().rstrip().split()
            clone_name = header[0]
            for line in inputs(f):
                barcode, read = line.split()[:2]
                read = int(read)
                if read > 0:
                    B2N[barcode] += 1

    N = len(datafiles) + 1
    valid_barcodes = sorted([barcode for barcode, n in B2N.items() if n == N])
    print(*valid_barcodes, sep='\n')

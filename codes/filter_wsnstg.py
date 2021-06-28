#!/usr/bin/env python3
import sys
import os
import re
from collections import defaultdict
from barcodelib import inputs

in_suffix = 'wsnstg'
out_suffix = 'filtered.wsnstg'
outdir = './'

if __name__ == '__main__':
    infiles = sys.argv[1:]
    population_barcodes = {barcode for barcode in inputs()}
    for infile_name in infiles:
        outfile_name = os.path.join(outdir,
            re.sub('\.'+in_suffix+'$', '', os.path.basename(infile_name)) + '.' + out_suffix
        )
        if os.path.exists(outfile_name):
            continue
        with open(infile_name) as infile:
            with open(outfile_name, 'w') as outfile:
                header = infile.readline().rstrip('\n')
                print(header, file=outfile)
                for line in inputs(infile):
                    wsn, stg, read = line.split()
                    if wsn in population_barcodes:
                        print(line, file=outfile)

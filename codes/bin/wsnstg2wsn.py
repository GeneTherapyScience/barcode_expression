#!/usr/bin/env python3
import sys
import os
from collections import defaultdict
import re
sys.path.append('../lib')
from alignment import seq_alignment
from barcodelib import inputs

usage = """\
Usage : $ python wsnstg2wsn.py [--outdir=<outdir>] [<file1> <file2> ...]\
"""

stg_template = "GGTGGCTTTACCAACAGTAC"
in_suffix = 'wsnstg'
out_suffix = 'wsn'

if __name__ == '__main__':
    infile_namelist = sys.argv[1:]
    diroption = "--outdir="
    if infile_namelist and infile_namelist[0][:len(diroption)] == diroption:
        outdir = infile_namelist[0][len(diroption):]
        infile_namelist = infile_namelist[1:]
    else:
        outdir = "./"
    if not infile_namelist:
        infile_namelist.append(None)

    for infile_name in infile_namelist:
        if infile_name is None:
            infile = sys.stdin
        else:
            outfile_name = os.path.join(outdir,
                re.sub('\.'+in_suffix+'$', '', os.path.basename(infile_name)) + '.' + out_suffix
                )
            if os.path.exists(outfile_name):
                continue
            infile = open(infile_name)

        wsn_data = defaultdict(lambda: [0,0]) # read, var
        header = infile.readline().rstrip('\n')
        for line in inputs(infile):
            wsn, stg, r = line.split()
            r = int(r)
            data = wsn_data[wsn]
            data[0] += r
            if stg != stg_template:
                data[1] += r

        if infile != sys.stdin:
            infile.close()
            outfile = open(outfile_name, 'w')
        else:
            outfile = sys.stdout

        print('WSN', 'Total', 'Mut', 'Mut_total', sep='\t', file=outfile)
        for wsn, data in sorted(wsn_data.items(), key=lambda x: (-x[1][0], -x[1][1], x[0])):
            r, mut = data
            if r > 0:
                mut_total = mut/r
            else:
                mut_total = 'N/A'
            print(wsn, r, mut, mut_total, sep='\t', file=outfile)

        if outfile != sys.stdout:
            outfile.close()

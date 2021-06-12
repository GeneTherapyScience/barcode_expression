#!/usr/bin/env python3
import sys
import os
from collections import defaultdict
import re
from alignment import seq_alignment
from merge_readerror import inputs

usage = """\
Usage : $ python wsnstg2stg.py [--outdir=<outdir>] [<file1> <file2> ...]\
"""

stg_template = "GGTGGCTTTACCAACAGTAC"
in_suffix = 'wsnstg'
out_suffix = 'stgmixd'

test = False
if test:
    # target_mixd = 'MXMXXXMXMXIIMXMXMMMXMXIIIIIIIII'
    stg_target = 'GATTCATCTCATCTATCAGAAAATAAATAAA'
    mixd, _ = seq_alignment(stg_template, stg_target)
    exit()

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

    letters = ['z'] + [chr(ord('a')+i) for i in range(20)]
    for infile_name in infile_namelist:
        if infile_name is None:
            infile = sys.stdin
        else:
            infile = open(infile_name)

        stg_data = dict() # read, var, mixd, distance
        header = infile.readline().rstrip('\n')
        for line in inputs(infile):
            wsn, stg, r = line.split()
            r = int(r)
            if stg in stg_data:
                data = stg_data[stg]
                data[0] += r
                data[1] += 1
            else:
                stg_data[stg] = [r,1] + list(seq_alignment(stg_template, stg))

        if infile != sys.stdin:
            infile.close()
            outfile_name = os.path.join(outdir,
                re.sub('\.'+in_suffix+'$', '', os.path.basename(infile_name)) + '.' + out_suffix
                )
            outfile = open(outfile_name, 'w')
        else:
            outfile = sys.stdout

        print('distance', 'reads', 'WSN\'s', 'stg', 'mixd', file=outfile)
        for stg, data in sorted(stg_data.items(), key=lambda x: (x[1][3], -x[1][0], -x[1][1], x[1][2], x[0])):
            r, v, mixd, d = data
            print(d, r, v, stg, mixd, sep='\t', file=outfile)

        if outfile != sys.stdout:
            outfile.close()

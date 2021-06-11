#!/usr/bin/env python3
import sys
import os
import re
from collections import defaultdict
from judge_stg import read_args
from math import ceil

usage = """\
Usage : $ python wsnstg2avarage.py --dictionary=<dictionary> [--outdir=<outdir>] [<file1> <file2> ...]\
"""

in_suffix = 'wsnstg'
out_suffix = 'mutinfo'

threshold_share = 10**(-5)

def wsn_thres(A, ratio=threshold_share):
    A = sorted(A)
    S = sum(A)
    prev = 0
    for a in A:
        if a == prev:
            S -= a
        elif a < S*ratio:
            S -= a
            prev = a
        else:
            return a
    else:
        return a

def inputs(f=sys.stdin):
    while True:
        line = f.readline()
        if not line: # EOF
            break
        else:
            yield line

if __name__ == '__main__':
    args, options = read_args(sys.argv)
    option_defaults = [
        ('dictionary', None),
        ('outdir', './')
    ]
    for op, def_val in option_defaults:
        globals()[op] = options[op] if op in options else def_val

    if dictionary is None:
        print('Give --dictionary.', file=sys.stderr)
    distance_dict = dict()
    with open(dictionary) as f:
        for line in inputs(f):
            stg, mixd, d = line.split()
            d = float(d)
            distance_dict[stg] = d

    infile_namelist = args[1:]
    if not infile_namelist:
        infile_namelist.append(None)

    for infile_name in infile_namelist:
        if infile_name is None:
            infile = sys.stdin
            outfile = sys.stdout
        else:
            outfile_name = os.path.join(outdir,
                re.sub('\.'+in_suffix+'$', '', os.path.basename(infile_name)) + '.' + out_suffix
                )
            infile = open(infile_name)
            outfile = open(outfile_name, 'w')

        header = infile.readline().rstrip('\n')
        output = [
            [''], ['total'], ['mut'], ['mut%'],
            ['mean distance (all):'], ['mean distance (mut):']]

        data = [line.split() for line in inputs(infile)]
        wsn_read = defaultdict(int)
        for i in range(2):
            totalread = 0
            mutread = 0
            Sd = 0
            mutvar = 0
            wsn_mut = set()
            wsn_all = set()
            for record in data:
                wsn, stg, r = record[:3]
                r = int(r)
                if i == 0:
                    wsn_read[wsn] += r
                elif wsn not in surviving_wsn:
                    continue

                d = distance_dict[stg]
                totalread += r
                wsn_all.add(wsn)
                if d > 0:
                    mutread += r
                    if wsn not in wsn_mut:
                        wsn_mut.add(wsn)
                        mutvar += 1
                Sd += d*r            
            if i == 0:
                output[0].extend(['read', 'WSN\'s'])
            else:
                output[0].extend(['read (>E-5)', 'WSN\'s (>E-5)'])
            output[1].extend([totalread, len(wsn_all)])
            output[2].extend([mutread, len(wsn_mut)])
            output[3].extend([
                round(mutread*100/totalread, 3),
                round(len(wsn_mut)*100/len(wsn_all), 3)])
            output[4].extend([round(Sd/totalread, 3), ''])
            output[5].extend([round(Sd/mutread, 3), ''])

            if i == 0:
                th = wsn_thres(wsn_read.values())
                surviving_wsn = {wsn for wsn, r in wsn_read.items() if r>=th}
        for line in output:
            print(*line, sep='\t', file=outfile)

        infile.close()
        outfile.close()

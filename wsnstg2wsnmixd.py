#!/usr/bin/env python3
import sys
import os
from collections import defaultdict
from alignment import seq_alignment
from merge_readerror import inputs

usage = """\
Usage : $ python wsnstg2wsnmixd.py [--outdir=<outdir>] [<file1> <file2> ...]\
"""

stg_template = "GGTGGCTTTACCAACAGTAC"
in_suffix = '.wsnstg'
out_suffix = '.wsnmixd'

def analyze_mixd(mixd):
    Is = []
    Ds = []
    p = -1
    i = 0
    d = 0
    for s in mixd:
        if s == 'I':
            i += 1
            d = 0
        else:
            p += 1
            if i > 0:
                Is.append((p,i))
                i = 0
            if s == 'D':
                d += 1
            elif d > 0:
                Ds.append((p-d,d))
                d = 0
    return Is, Ds

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

    letters = [chr(ord('a')+i) for i in range(20)] + ['z']
    for infile_name in infile_namelist:
        if infile_name is None:
            infile = sys.stdin
        else:
            infile = open(infile_name)

        mixd_data = dict()
        stg2mixd = dict()
        reads = defaultdict(int)
        ins_num = [0]*21
        del_num = [0]*20
        ins_len = [0]*31
        del_len = [0]*31
        header = infile.readline().rstrip('\n')
        for line in inputs(infile):
            wsn, stg, r = line.split()
            r = int(r)
            if stg in stg2mixd:
                mixd, distance = stg2mixd[stg]
            else:
                mixd, distance = seq_alignment(stg_template, stg)
                stg2mixd[stg] = (mixd, distance)
            reads[(wsn, mixd)] += r
            if mixd not in mixd_data:
                mixd_data[mixd] = (distance,) + analyze_mixd(mixd)
            _, Is, Ds = mixd_data[mixd]
            for p, i in Is:
                ins_num[p] += i*r
                if i > 30:
                    i = 30
                ins_len[i] += r
            for p, d in Ds:
                for x in range(p, p+d):
                    try:
                        del_num[x] += r
                    except:
                        print('error')
                        print('mixd', mixd)
                        print(x,p,d)
                        exit()
                if d > 30:
                    d = 30
                del_len[d] += r

        if infile != sys.stdin:
            infile.close()
            outfile_name = os.path.join(outdir,
                os.path.basename(infile_name).rstrip(in_suffix) + out_suffix
                )
            outfile = open(outfile_name, 'w')
        else:
            outfile = sys.stdout

        print('# reads', 'wsn', 'mixd', 'distance', 'ins', 'del', file=outfile)
        for pair, r in sorted(reads.items(), key=lambda x: (-x[1], x[0][0], mixd_data[x[0][1]][0], x[0][1])):
            wsn, mixd = pair
            distance, Is, Ds = mixd_data[mixd]
            Istr = ''.join([letters[p]+str(i) for p, i in Is])
            Dstr = ''.join([letters[p]+str(d) for p, d in Ds])
            print(r, wsn, mixd, distance, Is, Ds, sep='\t', file=outfile)

        print(infile_name)
        print('ins at each position (0-20):', ins_num[-1:]+ins_num[:-1], sep='\t')
        print('del at each position (1-20):', del_num[-1:]+del_num[:-1], sep='\t')
        print('ins of each lengths: (1-30+)', ins_len[1:], sep='\t')
        print('del of each lengths: (1-30+)', del_len[1:], sep='\t')

        if outfile != sys.stdout:
            outfile.close()

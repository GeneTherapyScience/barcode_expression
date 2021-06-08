#!/usr/bin/env python3
import sys
import os
from collections import defaultdict
import re
from alignment import seq_alignment
from merge_readerror import inputs

usage = """\
Usage : $ python wsnstg2wsnmixd.py [--outdir=<outdir>] [<file1> <file2> ...]\
"""

stg_template = "GGTGGCTTTACCAACAGTAC"
in_suffix = 'wsnstg'
out_suffix = 'wsnmixd'

def analyze_mixd(mixd):
    Is = []
    Ds = []
    p = 0
    i = d = 0
    for s in mixd +'M':
        if i > 0 and s != 'I':
            Is.append((p,i))
            i = 0
        if d > 0 and s != 'D':
            Ds.append((p-d+1,d))
            d = 0
        if s == 'I':
            i += 1
        elif s == 'D':
            d += 1
            p += 1
        else:
            p += 1
    return Is, Ds

test = False
if test:
    # target_mixd = 'MXMXXXMXMXIIMXMXMMMXMXIIIIIIIII'
    stg_target = 'GATTCATCTCATCTATCAGAAAATAAATAAA'
    mixd, _ = seq_alignment(stg_template, stg_target)
    print(analyze_mixd(mixd))
    print(analyze_mixd('MMDDDDDDDDDDDDDDDDDD'))
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

        mixd_data = dict()
        stg2mixd = dict()
        reads = defaultdict(int)
        ins_num = [0]*21
        del_num = [0]*21
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
                    del_num[x] += r
                if d > 30:
                    d = 30
                del_len[d] += r

        if infile != sys.stdin:
            infile.close()
            outfile_name = os.path.join(outdir,
                re.sub('\.'+in_suffix+'$', '', os.path.basename(infile_name)) + '.' + out_suffix
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
            print(r, wsn, mixd, distance, Istr, Dstr, sep='\t', file=outfile)

        print(infile_name)
        print('ins at each position (0-20):', ins_num, sep='\t')
        print('del at each position (1-20):', del_num[1:], sep='\t')
        print('ins of each lengths (1-30+):', ins_len[1:], sep='\t')
        print('del of each lengths (1-30+):', del_len[1:], sep='\t')

        if outfile != sys.stdout:
            outfile.close()

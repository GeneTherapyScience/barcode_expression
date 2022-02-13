#!/usr/bin/env python3
import sys
import os
import re
from judge_stg import read_args

usage = """\
Usage : $ python wsnstg2white.py --whitelist=<whitelist> [--outdir=<outdir>] [<file1> <file2> ...]\
"""

in_suffix = 'wsnstg'
out_suffix = 'white40.wsnstg'

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
        ('whitelist', None),
        ('outdir', './')
    ]
    for op, def_val in option_defaults:
        globals()[op] = options[op] if op in options else def_val

    if whitelist is None:
        print('Give --whitelist.', file=sys.stderr)
    whitelist_set = set()
    with open(whitelist) as f:
        whitelist_set = {line.strip() for line in inputs(f)}

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
        print(header, file=outfile)
        for line in inputs(infile):
            line = line.rstrip('\n')
            wsn, stg, reads = line.split()[:3]
            if stg in whitelist_set:
                print(line, file=outfile)
        infile.close()
        outfile.close()

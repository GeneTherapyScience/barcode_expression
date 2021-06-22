#!/usr/bin/env python3
import sys
import os
from judge_stg import read_args

usage = """\
Usage : $ python pickwsnlen.py [--outdir=<outdir>] [<file1> <file2> ...]\
"""

length_min = 30
length_max = 30

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
        ('wsncol', 1, int),
        ('outdir', './', str),
        ('noheader', False, bool),
    ]
    for op, def_val, dtype in option_defaults:
        globals()[op] = dtype(options[op]) if op in options else def_val
    infile_namelist = args[1:]
    if not infile_namelist:
        infile_namelist.append(None)

    for infile_name in infile_namelist:
        if infile_name is None:
            infile = sys.stdin
            outfile = sys.stdout
        else:
            suffix = os.path.basename(infile_name).split('.')[-1]
            outfile_name = os.path.join(outdir,
                os.path.basename(infile_name).rstrip(suffix) + "{}{}.{}".format(length_min,length_max,suffix)
                )
            infile = open(infile_name)
            outfile = open(outfile_name, 'w')

        if not noheader:
            header = infile.readline().rstrip('\n')
            print(header, file=outfile)
        for line in inputs(infile):
            wsn = line.split()[wsncol]
            if length_min <= len(wsn) <= length_max:
                print(line.rstrip('\n'), file=outfile)
        infile.close()
        outfile.close()

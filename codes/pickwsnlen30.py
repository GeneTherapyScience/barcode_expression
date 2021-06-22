#!/usr/bin/env python3
import sys
import os

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
            outfile = sys.stdout
        else:
            suffix = os.path.basename(infile_name).split('.')[-1]
            outfile_name = os.path.join(outdir,
                os.path.basename(infile_name).rstrip(suffix) + "{}{}.{}".format(length_min,length_max,suffix)
                )
            infile = open(infile_name)
            outfile = open(outfile_name, 'w')

        # header = infile.readline().rstrip('\n')
        # print(header, file=outfile)
        for line in inputs(infile):
            n, wsn = line.split()[:2]
            if length_min <= len(wsn) <= length_max:
                print(line.rstrip('\n'), file=outfile)
        infile.close()
        outfile.close()

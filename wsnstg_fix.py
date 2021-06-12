#!/usr/bin/env python3
import sys
import os
import shutil

usage = "Usage : $ python wsnstg_fix.py <file1> [<file2> ...]"

def inputs(f=sys.stdin):
    while True:
        line = f.readline()
        if not line: # EOF
            break
        else:
            yield line

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(usage, file=sys.stderr)
        exit()
    else:
        filelist = sys.argv[1:]

    tmpfilename = "/tmp/{}.tmp.{}".format(os.path.basename(sys.argv[0]), os.getpid())

    for infilename in filelist:
        if not os.path.exists(infilename):
            print("Skipping the file '{}', which does not exist...".format(infilename))
        with open(infilename) as in_file:
            with open(tmpfilename, 'w') as out_file:
                for line in inputs(in_file):
                    record = line.split()
                    record[0] = record[0].split('_')[-1]
                    print(*record, sep='\t', file=out_file)
        shutil.move(tmpfilename, infilename)

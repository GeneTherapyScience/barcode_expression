#!/usr/bin/env python3
import sys
import os
from collections import defaultdict

usage = "Usage : $ python sam2wsnstg_batch.py <listfile> <outdir>"

def inputs(f=sys.stdin):
    while True:
        line = f.readline()
        if not line: # EOF
            break
        else:
            yield line

def sam2wsnstg(in_filename, out_filename):
    with open(in_filename) as in_f:
        pairs = defaultdict(int)
        for line in inputs(in_f):
            if line[0] == '@':
                continue
            record = line.split()
            barcode = record[0].split("_TAG_")[-1]
            stg = record[9]
            pairs[(barcode, stg)] += 1

    with open(out_filename, "w") as out_f:
        print("WSN", "stg", "reads", sep="\t", file=out_f)
        for k, v in sorted(pairs.items(), key=lambda x: (x[0][0], -x[1], x[0][1])):
            print(*k, v, sep='\t', file=out_f)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(usage, file=sys.stderr)
        exit()
    else:
        samfilelist_name = sys.argv[1]
        outdir = sys.argv[2]

    with open(samfilelist_name) as list_f:
        samfile_list = list_f.readlines()

    for samfile_name in samfile_list:
        samfile_name = samfile_name.split('#')[0].strip()
        if not samfile_name:
            continue
        else:
            print("Processing input-file '{}'....".format(samfile_name), file=sys.stderr)
            outfile_name = os.path.join(outdir, os.path.basename(samfile_name).rstrip(".sam") + ".wsnstg")
            sam2wsnstg(samfile_name, outfile_name)
            print("Output-file '{}' was created.".format(outfile_name), file=sys.stderr)
    print("Completed.", file=sys.stderr)

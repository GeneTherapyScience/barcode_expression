#!/usr/bin/env python3
###
# usage: # python3 judge_stg.py [--reference=REFERENCE_SEQUENCE] blacklist_file < input_sequences > results
###
import os
import sys
from tqdm import tqdm, trange
import re
from barcodelib import inputs
from alignment import seq_distance

def read_args(argv):
    start_option='--'
    split_option='='
    args = []
    op_regex = re.compile(start_option+'(.+)'+split_option+'(.+)')
    options = dict()
    for S in argv:
        match = op_regex.fullmatch(S)
        if match:
            options[match.group(1)] = match.group(2)
        else:
            args.append(S)
    return args, options

default_reference = "GGTGGCTTTACCAACAGTAC"

if __name__ == '__main__':
    args, options = read_args(sys.argv)
    option_defaults = [
        ('reference', default_reference)
    ]
    for op, def_val in option_defaults:
        globals()[op] = options[op] if op in options else def_val
    blacklist_file = args[1]
    with open(blacklist_file) as f:
        blacklist_seq = [seq.strip() for seq in f.readlines()]
    input_sequences = [seq.strip() for seq in inputs()]
    N = len(input_sequences)
    print('stg', 'result', sep='\t')
    for i in trange(N):
        seq = input_sequences[i]
        ref_d = seq_distance(reference, seq)
        for x in blacklist_seq:
            if seq_distance(x, seq) < ref_d:
                print(seq, 'NG', sep='\t')
                break
        else:
            print(seq, 'OK', sep='\t')

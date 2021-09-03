from barcodelib import inputs
import numpy as np

def mutpos_mixd(mixd):
    ins_list, del_list, mis_list = [], [], []
    nI = 0
    for s in mixd:
        if s == 'I':
            nI += 1
        else:
            ins_list.append(nI)
            nI = 0
            if s == 'D':
                del_list.append(1)
                mis_list.append(0)
            elif s == 'X':
                del_list.append(0)
                mis_list.append(1)
            else:
                del_list.append(0)
                mis_list.append(0)
    else:
        ins_list.append(nI)

    return np.array(
        [1] + ins_list + del_list + mis_list
        + list(map(lambda x: int(bool(x)), ins_list)),
        dtype=int)

def stgmutpos(stgmixd_dict_file):
    mixdmutpos_data = dict()
    stgmutpos_data = dict()
    with open(stgmixd_dict_file) as f:
        for line in inputs(f):
            stg, mixd = line.split()[:2]
            if not mixd in mixdmutpos_data:
                mixdmutpos_data[mixd] = mutpos_mixd(mixd)
            stgmutpos_data[stg] = mixdmutpos_data[mixd]
    return stgmutpos_data

def pos_histogram(stgmutpos_data, wsnstg_filename):
    W = 1 + 21 + 20 + 20 + 21
    result = np.zeros(W, dtype=int)
    with open(wsnstg_filename) as f:
        f.readline()
        for line in inputs(f):
            wsn, stg, read = line.split()
            read = int(read)
            result += stgmutpos_data[stg] * read
    return result

import sys
import argparse
import re
if __name__ == '__main__':
    stgmixd_dict_file = "../stginfo/whitelist.sorted.stgmixd"

    parser = argparse.ArgumentParser()
    parser.add_argument('--perread', action='store_true',
                        help='output counts per read.')
    parser.add_argument('--excess', action='store_true',
                        help='output counts per read, over the plasmid sample.')
    parser.add_argument('--ratio', action='store_true',
                        help='output ratio to plasmid.')
    args = parser.parse_args()

    stgmutpos_dict = stgmutpos(stgmixd_dict_file)

    header = (
        ["I{}".format(i) for i in range(21)]
        + ["D{}".format(i) for i in range(1,21)]
        + ["X{}".format(i) for i in range(1,21)]
        + ["J{}".format(i) for i in range(21)]
    )

    if args.excess:
        filenames = list(map(lambda x: x.strip(), sys.stdin.readlines()))
        for filename in filenames:
            if re.match("Plasmid", filename):
                referencename = filename
                baseline = pos_histogram(stgmutpos_dict, filename)
                baseline = baseline[1:]/baseline[0]
                break
        else:
            print("Error: no reference file.", file=sys.stderr)
            exit()

        print("sample", *header, sep='\t')
        print(referencename, *baseline, sep='\t')
        for filename in filenames:
            if filename == referencename:
                continue
            hist = pos_histogram(stgmutpos_dict, filename)
            print(filename, *(hist[1:]/hist[0] - baseline), sep='\t')
    if args.ratio:
        filenames = list(map(lambda x: x.strip(), sys.stdin.readlines()))
        for filename in filenames:
            if re.match("Plasmid", filename):
                referencename = filename
                baseline = pos_histogram(stgmutpos_dict, filename)
                baseline = baseline[1:]/baseline[0]
                break
        else:
            print("Error: no reference file.", file=sys.stderr)
            exit()

        def ratio(data):
            N = len(data)
            ret = np.empty(N)
            for i in range(N):
                if baseline[i] == 0:
                    ret[i] = -1
                else:
                    ret[i] = data[i]/baseline[i]
            return np.round(ret, 3)

        print("sample", *header, sep='\t')
        print(referencename, *baseline, sep='\t')
        for filename in filenames:
            if filename == referencename:
                continue
            hist = pos_histogram(stgmutpos_dict, filename)
            print(filename, *ratio(hist[1:]/hist[0]), sep='\t')
    elif args.perread:
        print("sample", *header, sep='\t')
        for filename in inputs():
            hist = pos_histogram(stgmutpos_dict, filename)
            print(filename, *(hist[1:]/hist[0]), sep='\t')
    else:
        print("sample", "reads", *header, sep='\t')
        for filename in inputs():
            hist = pos_histogram(stgmutpos_dict, filename)
            print(filename, *hist, sep='\t')

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

    return np.array(ins_list + del_list + mis_list, dtype=int)

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
    W = 61
    result = np.zeros(W, dtype=int)
    with open(wsnstg_filename) as f:
        f.readline()
        for line in inputs(f):
            wsn, stg, read = line.split()
            read = int(read)
            result += stgmutpos_data[stg] * read
    return result

if __name__ == '__main__':
    stgmixd_dict_file = "../stginfo/whitelist.sorted.stgmixd"
    stgmutpos_dict = stgmutpos(stgmixd_dict_file)

    header = ["I{}".format(i) for i in range(21)] + ["D{}".format(i) for i in range(1,21)] + ["X{}".format(i) for i in range(1,21)]
    print("sample", *header, sep='\t')
    for filename in inputs():
        hist = pos_histogram(stgmutpos_dict, filename)
        print(filename, *hist, sep='\t')

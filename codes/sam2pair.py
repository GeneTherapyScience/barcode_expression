#!/user/bin/env python3
from merge_readerror import inputs
# from alignment import *
# from random import randrange
from collections import Counter, defaultdict
# from tqdm import tqdm, trange

stg_template = 'GGTGGCTTTACCAACAGTAC'

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # args = parser.parse_args()

    pairs = defaultdict(int)
    for line in inputs():
        if line[0] == '@':
            continue
        record = line.split()
        barcode = record[0].split('_TAG_')[-1]
        stg = record[9]
        pairs[(barcode, stg)] += 1

    print('WSN', 'stg', 'reads', sep='\t')
    for k, v in sorted(pairs.items(), key=lambda x: (x[0][0], -x[1], x[0][1])):
        print(*k, v, sep='\t')

#!/user/bin/env python3
from barcodelib import *
from alignment import *
from collections import Counter, defaultdict
from tqdm import tqdm, trange

stg_template = 'GGTGGCTTTACCAACAGTAC'

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # args = parser.parse_args()
    print(input())
    for line in inputs():

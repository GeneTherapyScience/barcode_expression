#!/usr/bin/env python3
from merge_readerror import *
from random import randrange
from collections import Counter
from tqdm import tqdm, trange

input()
total = single = 0
while True:
    try:
        line = input().strip()
        readnum = int(line.split()[1])
        total += readnum
        if readnum <= 1:
            single += readnum
    except EOFError:
        break

print('{}/{} = {}'.format(single, total, single/total))

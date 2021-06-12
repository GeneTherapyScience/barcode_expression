#!/usr/bin/env python3
import sys
from random import randrange
v = int(sys.argv[1])
n = int(sys.argv[2])
T = [0]*v
for _ in range(n):
    T[randrange(v)] += 1
print(sum(1 for t in T if t >= 2))

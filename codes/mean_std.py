import sys
import numpy as np
data = np.array([float(line) for line in sys.stdin.readlines()])
print(
    'mean = {}'.format(data.mean()),
    'std = {}'.format(data.std()),
    sep='\t')

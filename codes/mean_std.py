import sys
import numpy as np
data = np.array([float(line) for line in sys.stdin.readlines()])
print('mean:', data.mean(), sep='\t')
print('std:', data.std(), sep='\t')

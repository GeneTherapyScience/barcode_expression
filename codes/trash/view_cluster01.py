import sys
from barcodelib import inputs
import numpy as np

print('#', *sys.argv)
if len(sys.argv) < 2:
    inputfile = sys.stdin
else:
    inputfile = sys.argv[1]
    f = open(inputfile)

f.readline()
f.readline()

out_columns = ['Day0', 'Day7', 'Day14', 'NT', 'DTX', 'Sp', 'Day7NT', 'Day14NT', 'Day7DTX', 'Day14DTX', 'Day7Sp', 'Day14Sp']
members = [
    [1,1,1] + [0]*18,
    [0]*3 + [1]*9 + [0]*9,
    [0]*12 + [1]*9,
    [1]*3 + [1,1,1, 0,0,0, 0,0,0]*2,
    [0]*3 + [0,0,0, 1,1,1, 0,0,0]*2,
    [0]*3 + [0,0,0, 0,0,0, 1,1,1]*2,
    [0]*3 + [1]*3 + [0]*15,   # Day7NT
    [0]*12 + [1]*3 + [0]*6, # Day14NT
    [0]*6 + [1]*3 + [0]*12, # Day7DTX
    [0]*15 + [1]*3 + [0]*3, # Day14DTX
    [0]*9 + [1]*3 + [0]*9,  # Day7Sp
    [0]*18 + [1]*3,         # Day14Sp
]
M = 21

print(*out_columns, sep='\t')
output = []
for line in inputs(f):
    record = [
        round(sum(list(map(float,line.split()))[i]*members[w][i] for i in range(M)) / sum(members[w]),2)
            for w in range(len(out_columns))
    ]
    output.append(record)

output.sort(key=lambda x: x[4:], reverse=True)
for record in output:
    print(*record, sep='\t')

X = np.array(output)
print(np.corrcoef(X[:,4],X[:,5]))

#!/user/bin/python3
from merge_readerror import *

header, data = inputdata()
totalread = 0
totalvariation = 0
for record in data:
    barcode, read, altered = record[:3]
    totalread += read
    if read > 0:
        totalvariation += 1
print('variations', 'reads', sep='\t')
print(totalvariation, totalread, sep='\t')

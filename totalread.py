#!/user/bin/python3
from merge_readerror import *

header, data = inputdata()
totalread = 0
for record in data:
    barcode, read, altered = record[:3]
    totalread += read
print('variations', 'reads', sep='\t')
print(len(data), totalread, sep='\t')

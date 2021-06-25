from alignment import seq_alignment
from barcodelib import inputs
from tqdm import trange

default_reference = "GGTGGCTTTACCAACAGTAC"

lines = [line.strip() for line in inputs()]
N = len(lines)
for i in trange(N):
	stg = lines[i]
	mixd, distance = seq_alignment(default_reference, stg)
	print(stg, mixd, distance, sep='\t')

#!/usr/bin/env python3
import sys
sys.path.append('../lib')
from barcodelib import inputs

plasmid_file = "Plasmid_20nt345.merge.extracted.reformat.white40.ref.uniq"
barcode_cell_files = [
    "NANOG_112.merge.extracted.reformat.white40.ref.uniq",
    "NANOG_128.merge.extracted.reformat.white40.ref.uniq",
    "NANOG_151.merge.extracted.reformat.white40.ref.uniq",
    "POU5F1_112.merge.extracted.reformat.white40.ref.uniq",
    "POU5F1_114.merge.extracted.reformat.white40.ref.uniq",
    "SOX2_101.merge.extracted.reformat.white40.ref.uniq",
    "SOX2_98.merge.extracted.reformat.white40.ref.uniq",
]

output = "valid_barcodes.wsn"

plasmid_barcodes = set()
with open(plasmid_file) as f:
    f.readline()
    for line in inputs(f):
        _, barcode  = line.split()
        plasmid_barcodes.add(barcode)
cell_barcodes = set()
for barcode_file in barcode_cell_files:
    with open(barcode_file) as f:
        f.readline()
        for line in inputs(f):
            _, barcode = line.split()
            cell_barcodes.add(barcode)

cell_barcodes &= plasmid_barcodes
with open(output, 'w') as f:
    for barcode in sorted(cell_barcodes):
        print(barcode, file=f)

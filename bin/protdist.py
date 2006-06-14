#!/usr/bin/env python

from rasmus import phylip, fasta
import sys


files = sys.argv[1:]

for f in files:
    print f
    seqs = fasta.readFasta(f)
    phylip.protdist(seqs, f + ".dist")




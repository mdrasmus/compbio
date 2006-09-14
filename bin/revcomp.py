#!/usr/bin/env python


import sys
from rasmus import alignlib

seq = []
for line in sys.stdin:
    seq.append(line.rstrip())
seq = "".join(seq)

seq2 = alignlib.revcomp(seq)

for i in xrange(0, len(seq2), 60):
    print seq2[i:i+60]

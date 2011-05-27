#!/usr/bin/env python

import sys
from rasmus import fasta

fn = sys.argv[1]

names, seqs = zip(*fasta.read_fasta(fn))

for i in xrange(len(names)):
    if seqs[i][-1] == "*":
        seqs[i] = seqs[i][:-1]
    print ">" + names[i]
    print seqs[i]



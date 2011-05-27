#!/usr/bin/env python

import os, sys

from rasmus.bio import fasta


if len(sys.argv) < 3:
    print "fasta-get.py <fasta> <key>"
    sys.exit(1)

seqs = fasta.read_fasta(sys.argv[1])
print seqs[sys.argv[2]]




#!/usr/bin/env python

from rasmus import fasta, util
import sys


for fn in sys.argv[1:]:
    util.tic(fn)
    fa = fasta.read_fasta(fn, errors=False)

    fasta.write_fasta(fn, fa)
    util.toc()



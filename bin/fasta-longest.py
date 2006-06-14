#!/usr/bin/env python

from rasmus import fasta, util
import sys


for fn in sys.argv[1:]:
    util.tic(fn)
    fa = fasta.readFasta(fn, errors=False)

    fasta.writeFasta(fn, fa)
    util.toc()



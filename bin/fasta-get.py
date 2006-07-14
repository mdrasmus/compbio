#!/usr/bin/env python

import os, sys

from rasmus import fasta, util


if len(sys.argv) < 5:
    print "fasta-get.py <fasta> <key> <start> <end> [<strand>]"
    sys.exit(1)


if len(sys.argv) == 6:
    fafile, key, start, end, strand = sys.argv[1:6]
    strand = int(strand)
else:
    fafile, key, start, end = sys.argv[1:5]
    strand = 1

start = int(start)
end = int(end)


faindex = fasta.FastaIndex(fafile)
print faindex.get(key, start, end, strand)


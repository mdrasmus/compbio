#!/usr/bin/env python


import sys
from rasmus import genomeutil

seq = ""

for line in sys.stdin:
    seq += line.rstrip()

sys.stdout.write(genomeutil.reverseComplement(seq))

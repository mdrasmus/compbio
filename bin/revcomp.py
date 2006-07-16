#!/usr/bin/env python


import sys
from rasmus import alignlib

seq = ""

for line in sys.stdin:
    seq += line.rstrip()

sys.stdout.write(alignlib.revcomp(seq))

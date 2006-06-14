#!/usr/bin/env python
import sys
from rasmus import util


if (len(sys.argv) == 1):
    print "usage: cut.py [cols]"
    sys.exit(1)

# parse args
cols = map(int, sys.argv[1:])

for line in sys.stdin:
    tokens = line.rstrip().split("\t")
    print "\t".join(util.sublist(tokens, cols))




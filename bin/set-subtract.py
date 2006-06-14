#!/usr/bin/env python

from rasmus import util
import sys

if len(sys.argv) < 2:
    print >>sys.stderr, "usage: set-subtract.py <setA> <setB>"
    print >>sys.stderr, "    find  A - B"
    sys.exit(1)

set2 = util.makeset(util.readStrings(sys.argv[2]))

for line in file(sys.argv[1]):
    if line.rstrip() not in set2:
        print line,

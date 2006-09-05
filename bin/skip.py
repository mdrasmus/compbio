#!/usr/bin/env python

import sys

if (len(sys.argv) == 1):
    print "usage: skip.py [nlines]"
    sys.exit(1)

skip = int(sys.argv[1])


if skip > 0:
    n = 0
    for line in sys.stdin:
        n += 1
        if (n >= skip):
            break

for line in sys.stdin:
    print line,


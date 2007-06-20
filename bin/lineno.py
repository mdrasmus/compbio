#!/usr/bin/env python

import sys

if len(sys.argv) > 1:
    width = int(sys.argv[1])
    patt = "%%%dd  " % width
else:
    patt = "%d  "

i = 0
for line in sys.stdin:
    i += 1
    sys.stdout.write(patt % i)
    sys.stdout.write(line)
    sys.stdout.flush()

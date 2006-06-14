#!/usr/bin/env python

import sys


if len(sys.argv) < 4:
    print "usage: dnaget.py <start> <end> <file>"
    sys.exit(1)

start = int(sys.argv[1])
end = int(sys.argv[2])
filename = sys.argv[3]


infile = file(filename)
infile.seek(start - 1)

step = 10000
i = start

while i <= end:
    if i + step > end + 1:
        step = end + 1 - i
    i += step
    sys.stdout.write(infile.read(step))

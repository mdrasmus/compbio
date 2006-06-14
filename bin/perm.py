#!/usr/bin/python

import sys

if len(sys.argv) < 2:
    print "usage: perm.py <perm>"
    sys.exit(1)

infile = file(sys.argv[1])
perm = []
for line in infile:
    perm.append(int(line))
infile.close()

lines = []
for line in sys.stdin:
    lines.append(line)

for i in range(len(lines)):
    print lines[perm[i]],

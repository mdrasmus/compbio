#!/usr/bin/env python 

from rasmus import util
import sys

options = [
    ["l:", "line=", "line", "AUTO<line no>"],
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)

keys = map(int, param["line"][-1].split(","))
keyset = set(keys)
lines = {}

if len(rest) == 0:
    infile = sys.stdin
else:
    infile = file(rest[0])

found = 0
i = 0
for line in infile:
    if i in keyset:
        lines[i] = line
        if len(lines) == len(keys):
            break
    i += 1

for key in keys:
    print lines[key],

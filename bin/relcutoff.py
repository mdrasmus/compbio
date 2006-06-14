#!/usr/bin/python

import sys
import util

# parse args
if len(sys.argv) < 2:
    print "usage: relcutoff.py <cutoff>"
    sys.exit(1)

cutoff = float(sys.argv[1])


# read in data
mat = util.Dict(2)

for line in sys.stdin:
    (gene1, gene2, score) = line.split()[:3]    
    mat[gene1][gene2] = float(score)
    mat[gene2][gene1] = float(score)


# now print with cutoff
keys1 = mat.keys()
keys1.sort()
for key1 in keys1:
    row = mat[key1]
    keys = row.keys()
    top = max(row.values())
    
    keys = filter(lambda k: row[k] >= cutoff * top, keys)
    keys.sort()
    
    for key in keys:
        print key1, key, row[key]

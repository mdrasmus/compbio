#!/usr/bin/env python

import sys
from rasmus import fasta, util

options = [
    ["i:", "identity:", "identity", "AUTO<min identity stretch>"],
    ["n:", "num:", "num", "AUTO<min number of stretches>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)

identitylen = float(param["identity"][-1])
num = int(param["num"][-1])

for fn in rest:
    print "."
    fa = fasta.read_fasta(fn)
    conserve = fasta.calc_conservation(fa)
    
    matches = util.islands(conserve)
    
    if "*" not in matches:
        continue
    
    lens = [i[1] - i[0] for i in matches["*"]]
    
    lens.sort(lambda a,b: cmp(b, a))
    
    print lens[:20]

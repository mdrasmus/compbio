#!/usr/bin/python

import fasta
import util
import sys

fastafile, indexfile, labelfile = sys.argv[1:]


f = fasta.readFasta(fastafile)

keys = f.keys()
lookup = util.list2lookup(keys)

f2 = {}
for key in keys:
    f2[str(lookup[key])] = f[key]

util.writeVector(labelfile, keys)
fasta.writeFasta(indexfile, f2)

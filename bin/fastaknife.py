#!/usr/bin/env python

from rasmus import fasta, util
import sys


options = [
    ["k:", "keys=", "keys", "AUTO<keys>"],
    ["f:", "fasta=", "fasta", "AUTO<fasta file>"]
]



# check args
try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


keys = util.readStrings(param["keys"][-1])
fa  = fasta.readFasta(param["fasta"][-1])

fa2 = util.subdict(fa, keys)

fasta.writeFasta(sys.stdout, fa2)

#!/usr/bin/env python

from rasmus import phylip, fasta, util
import sys



options = [
  ["p:", "prog=", "prog", "<program name>",
    {"default": ["protdist"]}]
]


param = util.parseOptions(sys.argv, options, quit=True)



for f in param[""]:
    print f
    seqs = fasta.readFasta(f)
    
    if param["prog"][-1] == "protdist":
        phylip.protdist(seqs, f + ".dist")
    elif param["prog"][-1] == "dnadist":
        phylip.dnadist(seqs, f + ".dist")
    else:
        raise "unknown program '%s'" % param["prog"][-1]



    




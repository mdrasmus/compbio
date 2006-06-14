#!/usr/bin/env python

import sys
from rasmus import phylip, fasta, util, mrbayes


options = [
 ["f:", "fasta=", "fasta", "AUTO<fasta>"],
 ["l:", "label=", "label", "AUTO<output labels>"], 
 ["p:", "phylip=", "phylip", "AUTO<output phylip>"], 
 ["n:", "nexus=", "nexus", "AUTO<output nexus>"]
]


try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


if "phylip" in param and "label" in param:
    seqs = fasta.readFasta(param["fasta"][-1])
    labels = phylip.fasta2phylip(file(param["phylip"][-1], "w"), seqs)
    util.writeVector(param["label"][-1], labels)


if "nexus" in param:
    names, seqs = fasta.readFastaOrdered(param["fasta"][-1])
    mrbayes.writeNexus(file(param["nexus"][-1], "w"), names, seqs)
    


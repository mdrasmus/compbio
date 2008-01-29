#!/usr/bin/env python

import sys

from rasmus import util
from rasmus.bio import phylip, fasta, mrbayes


options = [
 ["f:", "fasta=", "fasta", "<fasta>"],
 ["l:", "label=", "label", "<output labels>"], 
 ["p:", "phylip=", "phylip", "<output phylip>"], 
 ["n:", "nexus=", "nexus", "<output nexus>"],
 ["", "nostrip", "nostrip", "",
    {"single": True}],
 ["t:", "seqtype=", "seqtype", "dna|pep",
  {"default": "dna",
   "single": True}],
]



param = util.parseOptions(sys.argv, options)

if "phylip" in param:
    seqs = fasta.readFasta(param["fasta"][-1])
    labels = phylip.fasta2phylip(file(param["phylip"][-1], "w"), seqs,
                                 stripNames=not param["nostrip"])
    
    if "label" in param:
        util.writeVector(param["label"][-1], labels)


if "nexus" in param:
    seqs = fasta.readFasta(param["fasta"][-1])
    mrbayes.writeNexus(file(param["nexus"][-1], "w"), seqs.keys(), seqs.values(),
                       format=param["seqtype"], seqwidth=util.INF)
    


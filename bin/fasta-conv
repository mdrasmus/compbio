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
    seqs = fasta.read_fasta(param["fasta"][-1])
    labels = phylip.write_phylip_align(
        file(param["phylip"][-1], "w"), seqs,
        strip_names=not param["nostrip"])
    
    if "label" in param:
        util.write_list(param["label"][-1], labels)


if "nexus" in param:
    seqs = fasta.read_fasta(param["fasta"][-1])
    mrbayes.write_nexus(file(param["nexus"][-1], "w"),
                        seqs.keys(), seqs.values(),
                        format=param["seqtype"])
    


#!/usr/bin/env python

from rasmus import util, phylip, muscle, clustalw, fasta

import sys

options = [
    ["b:", "bootiters=", "bootiters", "<# iterations>"],
]

param = util.parseOptions(sys.argv, options, quit=True)



if "bootiters" not in param:
    aln = fasta.readFasta(sys.stdin)
    tree = phylip.proml(aln)
    tree.writeNewick(sys.stdout)
else:
    aln = fasta.readFasta(sys.stdin)
    trees = phylip.bootProml(aln, int(param["bootiters"][-1]), jumble=1)
    
    tree = phylip.consense(trees)
    tree.writeNewick(sys.stdout)
    
    #for tree in trees:
    #    tree.writeNewick(sys.stdout)

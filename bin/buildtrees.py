#!/usr/bin/env python

from rasmus import util, phylip, muscle, clustalw, fasta, bionj

import sys


options = [
    ["b:", "boot=", "boot", "AUTO<iters>"],
    ["d:", "dir=", "dir", "AUTO<alignments directory>"],
    ["i:", "index=", "index", "AUTO<starting index>"]
    ]
    

try:
    param, rest = util.parseArgs(sys.argv, options, "alignments")
except:
    sys.exit(1)

if "boot" in param:
    iters = int(param["boot"][-1])

if "dir" in param:
    files = util.listFiles(param["dir"][-1], ".align")
else:
    files = rest

if "index" in param:
    index = int(param["index"][-1])
else:
    index = 0

for align in files[index:]:
    util.log(align)
    aln = fasta.readFasta(align)
    
    if len(aln) > 2:
        if "boot" in param:
            trees = phylip.bootNeighbor(aln, iters)
            phylip.writeBootTrees(align.replace(".align", ".boot.tree"), trees)
        else:
            tree = bionj.bionj(aln)
            tree.writeNewick(file(align.replace(".align", ".tree"), "w"))
    



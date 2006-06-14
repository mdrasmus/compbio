#!/usr/bin/env python

from rasmus import util, blast
import sys

options = [
    ["g:", "genes=", "genes", "AUTO<acceptable genes>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


if "genes" in param:
    genefilter = True
    genes = util.readStrings(param["genes"][-1])
    lookup = util.makeset(genes)
else:
    genefilter = False


best = util.Dict(1, [0,None])

for f in rest:
    reader = blast.BlastReader(f)
    for line in reader:
        name1 = blast.query(line)
        name2 = blast.subject(line)
        score = blast.bitscore(line)

        # only use acceptable genes
        if genefilter:
            if name1 not in lookup or name2 not in lookup:
                continue

        if score > best[name1][0]:
            best[name1] = [score, name2]
        if score > best[name2][0]:
            best[name2] = [score, name1]

for gene1 in best:
    gene2 = best[gene1][1]
    
    # see if the best is bidirectional
    # but only print one direction
    if gene1 == best[gene2][1] and gene1 < gene2:
        print "%s\t%s\t%f" % (gene1, gene2, best[gene1][0])
        

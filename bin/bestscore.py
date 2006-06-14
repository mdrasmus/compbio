#!/usr/bin/python

from rasmus import util, graph
import getopt
import sys



if len(sys.argv) < 2:
    print >>sys.stderr, "usage: bestscore.py [-m matches]"
    sys.exit(1)


params, rest = util.getopt(sys.argv[1:], "m:")

best = {}

# pull in matches above threshold
for matchfile in params["-m"]:
    print >>sys.stderr, matchfile
    for line in file(matchfile):
        (gene1, gene2, score) = line.split()[:3]
        score = float(score)
        
        if gene1 == gene2:
            continue
        
        if gene1 not in best or score > best[gene1][1]:
            best[gene1] = (gene2, score)
        if gene2 not in best or score > best[gene2][1]:
            best[gene2] = (gene1, score)

for gene1, (gene2, score) in best.items():
    print gene1, gene2, score

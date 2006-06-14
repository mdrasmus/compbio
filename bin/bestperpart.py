#!/usr/bin/env python

from rasmus import blast, util
import sys

options = [
    ["p:", "part=", "part", "AUTO<part file>"],
    ["c:", "cols=", "cols", "AUTO<fields>"]
]


try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


parts = util.readDelim(param["part"][-1])
gene1col, gene2col, scorecol = map(int, param["cols"][-1].split(","))

lookup = {}

for i in xrange(len(parts)):
    for gene in parts[i]:
        lookup[gene] = i

hits = util.Dict(2, [0, []])


for filename in rest:
    for line in util.DelimReader(filename):
        gene1 = line[gene1col]
        gene2 = line[gene2col]
        score = float(line[scorecol])
        
        if gene2 in lookup:
            part2 = lookup[gene2]        
            if score > hits[gene1][part2][0]:
                hits[gene1][part2] = (score, line)
        
        if gene1 in lookup:
            part1 = lookup[gene1]
            if score > hits[gene2][part1][0]:
                hits[gene2][part1] = (score, line)


for gene in hits:
    for part in hits[gene].values():
        print "\t".join(part[1])


